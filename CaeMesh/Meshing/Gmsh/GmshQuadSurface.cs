// PrePoMax - Copyright (C) 2016-2026 Matej Borovinšek
//
// Licensed under the terms defined in the LICENSE file located in the root directory of this source code.
//
// Source code: https://github.com/PrePoMax/PrePoMax
//
// Author: Matej Borovinšek
// Contributors:

using System;
using System.Collections.Generic;
using System.Linq;
using CaeGlobals;
using GmshCommon;

namespace CaeMesh
{
    /// <summary>
    /// A surface of the gmsh model as the surface a mesh is smoothed on in space: a node is moved onto it by the closest
    /// point of the surface and takes the normal of the surface there, and its parameters when those of the nodes
    /// were given. A closest point farther from the node than the node moved is none, as a move in the tangent plane
    /// leaves the surface by far less than its own length, and neither is one the parameters of the surface do not
    /// give back to a millionth of the size of the mesh, as a point beside the face.
    /// </summary>
    internal class GmshQuadSurface : IQuadSurface
    {
        // Variables                                                                                                    
        private const int NumOfGridSmoothingSweeps = 50;
        private const double ParameterTolerance = 1E-6;
        private readonly int _surfaceId;
        private readonly double[][] _uv;
        private readonly double _size;


        // Properties                                                                                                   
        public double[][] Uv { get { return _uv; } }


        // Constructors                                                                                                 
        public GmshQuadSurface(int surfaceId, double[][] uv, double size)
        {
            _surfaceId = surfaceId;
            _uv = uv;
            _size = size;
        }


        // Methods                                                                                                      
        /// <summary>
        /// Smooths the free nodes of the quads and triangles of a surface in space (QuadSmoothing.SmoothRegionInSpace),
        /// every other node staying where it is, and gives them their new coordinates and parameters. A corner at a
        /// pole node, where the surface has no normal, is measured against the normal of its own element. The
        /// smoothed mesh is kept only when its worst scaled Jacobian in space is not lower than before; otherwise
        /// nothing changes. Returns whether it was kept.
        /// </summary>
        public static bool Smooth(int surfaceId, double[][] xyz, double[][] uv, int[][] elements, int[] freeNodes,
                                  ICollection<int> poleNodes, int maxNumOfSweeps, double target)
        {
            if (freeNodes.Length == 0 || elements.Length == 0) return false;
            double orientation;
            double[][] normals = GetUnitNormals(surfaceId, xyz, uv, elements, poleNodes, out orientation);
            double before = QuadSmoothing.GetSpaceMinScaledJacobian(xyz, normals, elements, orientation);
            double[][] smoothedXyz = (double[][])xyz.Clone();
            GmshQuadSurface surface = new GmshQuadSurface(surfaceId, (double[][])uv.Clone(), GetSize(xyz));
            double after = QuadSmoothing.SmoothRegionInSpace(smoothedXyz, normals, elements, freeNodes, orientation,
                                                             maxNumOfSweeps, target, surface);
            if (!(after >= before)) return false;
            foreach (var node in freeNodes)
            {
                xyz[node] = smoothedXyz[node];
                uv[node] = surface.Uv[node];
            }
            return true;
        }
        /// <summary>
        /// The unit normals of the surface at the nodes, null at a pole node and where the surface has none, and the
        /// way the elements turn against them.
        /// </summary>
        internal static double[][] GetUnitNormals(int surfaceId, double[][] xyz, double[][] uv, int[][] elements,
                                                 ICollection<int> poleNodes, out double orientation)
        {
            double[][] normals = GmshQuadMesher.GetNormals(surfaceId, uv);
            for (int i = 0; i < normals.Length; i++)
            {
                double length = Math.Sqrt(normals[i][0] * normals[i][0] + normals[i][1] * normals[i][1] +
                                          normals[i][2] * normals[i][2]);
                if ((poleNodes != null && poleNodes.Contains(i)) || !(length > 0)) normals[i] = null;
                else for (int d = 0; d < 3; d++) normals[i][d] /= length;
            }
            double energy;
            double sum = 0;
            foreach (var element in elements)
            {
                for (int c = 0; c < element.Length; c++)
                {
                    if (normals[element[c]] != null)
                        sum += QuadSmoothing.GetSpaceCornerQuality(xyz, normals, element, c, 1, out energy);
                }
            }
            orientation = sum < 0 ? -1 : 1;
            return normals;
        }
        /// <summary>
        /// The worst scaled Jacobian in space of the corners of the quads and triangles of a surface, at a pole node
        /// against the normal of its own element, and the number of elements with a corner that turns the other way
        /// or has no area.
        /// </summary>
        public static double Judge(int surfaceId, double[][] xyz, double[][] uv, int[][] elements,
                                   ICollection<int> poleNodes, out int numOfFoldedElements)
        {
            double orientation;
            double[][] normals = GetUnitNormals(surfaceId, xyz, uv, elements, poleNodes, out orientation);
            double energy;
            double min = 1;
            numOfFoldedElements = 0;
            foreach (var element in elements)
            {
                double elementMin = 1;
                for (int c = 0; c < element.Length; c++)
                    elementMin = Math.Min(elementMin, QuadSmoothing.GetSpaceCornerQuality(xyz, normals, element, c,
                                                                                         orientation, out energy));
                if (!(elementMin > 0)) numOfFoldedElements++;
                min = Math.Min(min, elementMin);
            }
            return min;
        }
        /// <summary>
        /// Smooths the inner nodes of a grid a mesher laid and smoothed in a plane frame of the parameters of the
        /// surface, whose nodes are the boundary nodes followed by the inner ones, in space: the frame keeps the
        /// lengths and angles of a curved surface only in the mean, so the grid that is best in it is not the best on
        /// the surface. The frame of a plane is the plane itself, where nothing is left to do. Returns whether the
        /// inner nodes were moved.
        /// </summary>
        public static bool SmoothGrid(int surfaceId, double[][] boundaryXyz, double[][] boundaryUv, double[][] innerXyz,
                                      double[][] innerUv, int[][] elements, ICollection<int> poleNodes = null)
        {
            if (innerUv.Length == 0 || Gmsh.Model.GetType(2, surfaceId) == "Plane") return false;
            int numOfBoundaryNodes = boundaryXyz.Length;
            double[][] xyz = boundaryXyz.Concat(innerXyz).ToArray();
            double[][] uv = boundaryUv.Concat(innerUv).ToArray();
            int[] freeNodes = Enumerable.Range(numOfBoundaryNodes, innerUv.Length).ToArray();
            if (!Smooth(surfaceId, xyz, uv, elements, freeNodes, poleNodes, NumOfGridSmoothingSweeps,
                        SubmapFaceDetector.MinScaledJacobian)) return false;
            for (int n = 0; n < innerUv.Length; n++)
            {
                innerXyz[n] = xyz[numOfBoundaryNodes + n];
                innerUv[n] = uv[numOfBoundaryNodes + n];
            }
            return true;
        }
        /// <summary>
        /// SmoothGrid for a grid that is yet to be judged: on a curved surface the grid is smoothed in space whatever
        /// the frame made of it, and the number of its folded quads and its worst scaled Jacobian become those it has
        /// in space (Judge), since a grid the frame folds or squeezes may lie well on the surface and one it likes may
        /// not. On a plane both stay as the frame gave them. Returns whether the grid was judged in space.
        /// </summary>
        public static bool SmoothGrid(int surfaceId, double[][] boundaryXyz, double[][] boundaryUv, double[][] innerXyz,
                                      double[][] innerUv, int[][] elements, ICollection<int> poleNodes,
                                      ref int numOfFoldedElements, ref double minScaledJacobian)
        {
            if (Gmsh.Model.GetType(2, surfaceId) == "Plane") return false;
            SmoothGrid(surfaceId, boundaryXyz, boundaryUv, innerXyz, innerUv, elements, poleNodes);
            minScaledJacobian = Judge(surfaceId, boundaryXyz.Concat(innerXyz).ToArray(),
                                      boundaryUv.Concat(innerUv).ToArray(), elements, poleNodes,
                                      out numOfFoldedElements);
            return true;
        }
        /// <summary>
        /// The length of the diagonal of the box around the nodes.
        /// </summary>
        public static double GetSize(double[][] xyz)
        {
            double[] low = new double[] { double.MaxValue, double.MaxValue, double.MaxValue };
            double[] high = new double[] { double.MinValue, double.MinValue, double.MinValue };
            foreach (var point in xyz)
            {
                for (int d = 0; d < 3; d++)
                {
                    low[d] = Math.Min(low[d], point[d]);
                    high[d] = Math.Max(high[d], point[d]);
                }
            }
            double size2 = 0;
            for (int d = 0; d < 3; d++) size2 += (high[d] - low[d]) * (high[d] - low[d]);
            return xyz.Length > 0 ? Math.Sqrt(size2) : 0;
        }
        public void Project(double[][] xyz, double[][] normals, IList<int> nodes, double[][] previousXyz)
        {
            if (nodes.Count == 0) return;
            double[] coor = new double[3 * nodes.Count];
            for (int i = 0; i < nodes.Count; i++) Array.Copy(xyz[nodes[i]], 0, coor, 3 * i, 3);
            double[] closest;
            double[] parameters = null;
            double[] flatNormals = null;
            double[] values = null;
            // A point beside the face has no parameters on it, of which gmsh warns point by point
            double verbosity = Gmsh.Option.GetNumber("General.Verbosity");
            Gmsh.Option.SetNumber("General.Verbosity", 1);
            bool valid;
            try
            {
                Gmsh.Model.GetClosestPoint(2, _surfaceId, coor, out closest, out _);
                valid = closest != null && closest.Length == coor.Length;
                if (valid)
                {
                    Gmsh.Model.GetParametrization(2, _surfaceId, closest, out parameters);
                    valid = parameters != null && parameters.Length == 2 * nodes.Count;
                }
                if (valid)
                {
                    Gmsh.Model.GetNormal(_surfaceId, parameters, out flatNormals);
                    Gmsh.Model.GetValue(2, _surfaceId, parameters, out values);
                    valid = flatNormals != null && flatNormals.Length == coor.Length && values != null &&
                            values.Length == coor.Length;
                }
            }
            finally
            {
                Gmsh.Option.SetNumber("General.Verbosity", verbosity);
            }
            double tolerance2 = ParameterTolerance * _size * ParameterTolerance * _size;
            for (int i = 0; i < nodes.Count; i++)
            {
                int node = nodes[i];
                bool taken = valid;
                double distance2 = 0;
                double move2 = 0;
                double normal2 = 0;
                double error2 = 0;
                for (int d = 0; d < 3 && taken; d++)
                {
                    double c = closest[3 * i + d];
                    double n = flatNormals[3 * i + d];
                    if (double.IsNaN(c) || double.IsInfinity(c) || double.IsNaN(n)) taken = false;
                    distance2 += (c - xyz[node][d]) * (c - xyz[node][d]);
                    if (previousXyz != null)
                        move2 += (xyz[node][d] - previousXyz[i][d]) * (xyz[node][d] - previousXyz[i][d]);
                    normal2 += n * n;
                    error2 += (values[3 * i + d] - c) * (values[3 * i + d] - c);
                }
                if (!taken || (previousXyz != null && distance2 > move2) || !(normal2 > 0) || !(error2 <= tolerance2))
                {
                    if (previousXyz != null) xyz[node] = previousXyz[i];
                    continue;
                }
                double length = Math.Sqrt(normal2);
                xyz[node] = new double[] { closest[3 * i], closest[3 * i + 1], closest[3 * i + 2] };
                normals[node] = new double[] { flatNormals[3 * i] / length, flatNormals[3 * i + 1] / length,
                                               flatNormals[3 * i + 2] / length };
                if (_uv != null) _uv[node] = new double[] { parameters[2 * i], parameters[2 * i + 1] };
            }
        }
    }
}
