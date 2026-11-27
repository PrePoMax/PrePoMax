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
using System.Text.RegularExpressions;
using CaeGlobals;
using GmshCommon;

namespace CaeMesh
{
    /// <summary>
    /// Meshes a submapped face by SubmapMapping, the structured mapping every mesher shares, the way
    /// GmshThreeBlockMesher meshes a three-sided face: gmsh supplies the nodes its edges carry, their parameters on
    /// the surface and the evaluation of the surface, and the elements gmsh meshed on the face are replaced. A grid
    /// that cannot be made, folds or has a quad corner below SubmapFaceDetector.MinScaledJacobian leaves the face as
    /// gmsh meshed it; the grid of a curved surface, its rings included, is first smoothed on the surface itself and
    /// judged there, in space (GmshQuadSurface.SmoothGrid). On request the grid of a plane face without a ring
    /// keeps its lattice lines straight where that gives a good grid (SubmapMapping), and at right angles where it
    /// has no round hole and its boundary nodes lie on such lines. An edge that turns in its middle is split at its
    /// middle
    /// node into the two logical edges of the face, so it needs an even number of elements; one that turns twice at its nodes at a quarter and three quarters
    /// into three, so it needs a multiple of four. The nodes of the holes of a face are read loop by loop as those of
    /// its outer loop, and the holes are placed in the lattice by SubmapHolePlacement on the parameters of the nodes;
    /// a hole for which no place is found leaves the face as gmsh meshed it. The nodes of a round hole are read around
    /// it and split into its four logical edges from the node SubmapHolePlacement.GetRoundLoopOrder finds, so its edges
    /// must carry as many elements each and a multiple of four around it.
    /// A round hole with room beside it is meshed with a ring of elements around it: the grid is laid around a
    /// circle that lies a ring thickness outside the hole, node for node, and the ring joins the nodes of the hole to
    /// those of the circle. The nodes of the hole are then all shared by two elements, and the four corners of the
    /// square of the hole in the lattice, where five elements meet, lie on the circle and not on the hole. The ring is
    /// as thick as its elements are long along the hole, and no thicker than half of what lies between the hole and
    /// the nearest other loop of the face, a third towards a round hole, which may carry a ring of its own; where
    /// less than MinRingRoom of the element length lies there the hole is meshed without a ring.
    /// </summary>
    public static class GmshSubmapMesher
    {
        // Variables
        /// <summary>
        /// The least distance between a round hole and the other loops of its face, as a share of the length of the
        /// elements along the hole, with which the hole is meshed with a ring.
        /// </summary>
        public const double MinRingRoom = 1.5;
        private static readonly Regex LoopEdge = new Regex(@"(?<![A-Za-z])[Ee]dge(s?) (\d+)(?: and (\d+))?",
                                                           RegexOptions.Compiled);


        // Methods                                                                                                                  
        /// <summary>
        /// The face carries gmsh ids: its surface, its vertices and its edges, those of its holes as well. Without quads,
        /// as for a surface that is not recombined, every quad is split into two triangles. False, with the reason, when
        /// the face is left as gmsh meshed it.
        /// </summary>
        public static bool Mesh(SubmapFace face, bool quads, out string reason)
        {
            return Mesh(face, quads, out reason, false);
        }
        internal static bool Mesh(SubmapFace face, bool quads, out string reason, bool keepOrthogonal)
        {
            SubmapGrid grid;
            List<IntPtr>[][] loopTags;
            double[][][][] loopXyz;
            double[][][][] loopUv;
            int[][] gridQuads;
            double[][] innerUv;
            double[][] innerXyz;
            if (!Map(face, keepOrthogonal, out grid, out loopTags, out loopXyz, out loopUv, out gridQuads, out innerUv,
                     out innerXyz, out reason)) return false;
            // The tags of the boundary nodes of all loops in loop order, and the parameters of the outer loop, which
            // tell the way the grid turns
            IntPtr[] boundaryNodeTags = new IntPtr[grid.NumberOfBoundaryNodes];
            double[][] outerUv = new double[grid.NumberOfOuterBoundaryNodes][];
            int edge = 0;
            for (int l = 0; l < loopTags.Length; l++)
            {
                for (int k = 0; k < loopTags[l].Length; k++, edge++)
                {
                    for (int i = 0; i < grid.EdgeCounts[edge]; i++)
                    {
                        int nodeId = grid.GetBoundaryNodeId(edge, i);
                        boundaryNodeTags[nodeId] = loopTags[l][k][i];
                        if (l == 0) outerUv[nodeId] = loopUv[l][k][i];
                    }
                }
            }
            GmshThreeBlockMesher.ReplaceSurfaceMesh(face.FaceId, boundaryNodeTags, outerUv, gridQuads, innerUv,
                                                    innerXyz, quads);
            return true;
        }
        /// <summary>
        /// Whether the grid of the face, which carries gmsh ids and whose edges gmsh has meshed, would be kept, and its
        /// worst scaled Jacobian of a quad corner in space against the normal of the surface, at a pole against the
        /// normal of its own quad, as GmshLensMesher judges a lens grid. Nothing is written into the surface.
        /// </summary>
        public static bool Evaluate(SubmapFace face, out double minScaledJacobian, out string reason)
        {
            return Evaluate(face, out minScaledJacobian, out reason, false);
        }
        internal static bool Evaluate(SubmapFace face, out double minScaledJacobian, out string reason,
                                      bool keepOrthogonal)
        {
            minScaledJacobian = double.NaN;
            SubmapGrid grid;
            List<IntPtr>[][] loopTags;
            double[][][][] loopXyz;
            double[][][][] loopUv;
            int[][] gridQuads;
            double[][] innerUv;
            double[][] innerXyz;
            if (!Map(face, keepOrthogonal, out grid, out loopTags, out loopXyz, out loopUv, out gridQuads, out innerUv,
                     out innerXyz, out reason)) return false;
            // The inner nodes are keyed by negative tags, which no node of gmsh carries
            IntPtr[] nodeTags = new IntPtr[grid.NumberOfBoundaryNodes + innerUv.Length];
            Dictionary<IntPtr, double[]> tagXyz = new Dictionary<IntPtr, double[]>();
            Dictionary<IntPtr, double[]> tagUv = new Dictionary<IntPtr, double[]>();
            List<double[]> cornerXyz = new List<double[]>();
            int edge = 0;
            for (int l = 0; l < loopTags.Length; l++)
            {
                for (int k = 0; k < loopTags[l].Length; k++, edge++)
                {
                    cornerXyz.Add(loopXyz[l][k][0]);
                    for (int i = 0; i < grid.EdgeCounts[edge]; i++)
                    {
                        int nodeId = grid.GetBoundaryNodeId(edge, i);
                        nodeTags[nodeId] = loopTags[l][k][i];
                        tagXyz[nodeTags[nodeId]] = loopXyz[l][k][i];
                        tagUv[nodeTags[nodeId]] = loopUv[l][k][i];
                    }
                }
            }
            for (int n = 0; n < innerUv.Length; n++)
            {
                IntPtr tag = new IntPtr(-1 - n);
                nodeTags[grid.NumberOfBoundaryNodes + n] = tag;
                tagXyz[tag] = innerXyz[n];
                tagUv[tag] = innerUv[n];
            }
            IntPtr[] quadNodeTags = gridQuads.SelectMany(quad => quad.Select(nodeId => nodeTags[nodeId])).ToArray();
            int numOfFoldedQuads;
            minScaledJacobian = GmshSeamGridMesher.GetMinScaledJacobian(face.FaceId, quadNodeTags, tagXyz, tagUv,
                                                                        GetPoleTags(face.FaceId, loopTags, cornerXyz),
                                                                        out numOfFoldedQuads);
            if (numOfFoldedQuads > 0)
            {
                reason = numOfFoldedQuads + " of the " + gridQuads.Length + " quads of its grid fold in space";
                return false;
            }
            return true;
        }
        /// <summary>
        /// The tags of the nodes at the start of the logical edges, cornerXyz, that are poles of the surface.
        /// </summary>
        private static HashSet<IntPtr> GetPoleTags(int surfaceId, List<IntPtr>[][] loopTags, List<double[]> cornerXyz)
        {
            HashSet<IntPtr> poleTags = new HashSet<IntPtr>();
            double[] low = new double[] { double.MaxValue, double.MaxValue, double.MaxValue };
            double[] high = new double[] { double.MinValue, double.MinValue, double.MinValue };
            foreach (var point in cornerXyz)
            {
                for (int c = 0; c < 3; c++)
                {
                    low[c] = Math.Min(low[c], point[c]);
                    high[c] = Math.Max(high[c], point[c]);
                }
            }
            double size = GmshPoleChart.GetDistance(low, high);
            int p;
            double poleV;
            double farV;
            int corner = 0;
            for (int l = 0; l < loopTags.Length; l++)
            {
                for (int k = 0; k < loopTags[l].Length; k++, corner++)
                {
                    if (GmshPoleChart.FindPole(surfaceId, cornerXyz[corner], size, out p, out poleV, out farV))
                        poleTags.Add(loopTags[l][k][0]);
                }
            }
            return poleTags;
        }
        /// <summary>
        /// The grid of the face on the nodes gmsh meshed on its edges, with the tags, coordinates and parameters of the
        /// nodes of its logical edges loop by loop. False, with the reason, when the grid cannot be made, folds or has
        /// a quad corner below SubmapFaceDetector.MinScaledJacobian. The quads and the inner nodes are those of the
        /// grid, and with a ring around a round hole those of the ring as well: the boundary nodes of the grid along
        /// the hole are the nodes of the hole, and the nodes of the circle the grid is laid around are inner nodes that
        /// follow those of the grid.
        /// </summary>
        private static bool Map(SubmapFace face, bool keepOrthogonal, out SubmapGrid grid,
                                out List<IntPtr>[][] loopTags, out double[][][][] loopXyz,
                                out double[][][][] loopUv, out int[][] quads, out double[][] innerUv,
                                out double[][] innerXyz, out string reason)
        {
            grid = null;
            quads = null;
            innerUv = null;
            innerXyz = null;
            if (!TryPlaceHoles(face, out loopTags, out loopXyz, out loopUv, out _, out int[][] holeOffsets, out reason))
                return false;
            int surfaceId = face.FaceId;
            SubmapFace[] loops = face.Loops;
            int[][] loopDirections = face.LoopLogicalEdgeDirections;
            bool curved = Gmsh.Model.GetType(2, surfaceId) != "Plane";
            bool keepStraight = keepOrthogonal && !curved;
            keepOrthogonal = keepStraight && loops.All(loop => !loop.IsRound);
            double[][][][] gridXyz;
            double[][][][] gridUv;
            GetRingLoops(surfaceId, loops, loopXyz, loopUv, out gridXyz, out gridUv);
            // The nodes a ring is laid on are moved with those of the grid, which then keeps no straight lines
            for (int l = 0; l < loops.Length; l++) keepStraight = keepStraight && ReferenceEquals(gridUv[l], loopUv[l]);
            try
            {
                grid = SubmapMapping.Map(loopDirections, gridUv, gridXyz, holeOffsets, uv =>
                {
                    double[] coor;
                    Gmsh.Model.GetValue(2, surfaceId, uv, out coor);
                    return coor;
                }, keepOrthogonal, keepStraight);
            }
            catch (CaeException ex)
            {
                reason = NameLoopEdges(face, ex.Message);
                return false;
            }
            // With a ring the grid is judged as it is meshed, its nodes moved together with those of the ring
            int numOfFoldedQuads = grid.NumberOfFoldedQuads;
            double minScaledJacobian = grid.MinScaledJacobian;
            AddRings(surfaceId, grid, loopXyz, loopUv, gridXyz, gridUv, out quads, out innerUv, out innerXyz,
                     ref numOfFoldedQuads, ref minScaledJacobian);
            // The grid of a curved surface, with its rings, is smoothed and judged on the surface itself; a pole of the
            // surface at a corner has no normal there
            if (curved)
            {
                HashSet<IntPtr> poleTags = GetPoleTags(surfaceId, loopTags, loopXyz.SelectMany(loop => loop)
                                                                                   .Select(side => side[0]).ToList());
                HashSet<int> poleNodes = new HashSet<int>();
                double[][] boundaryXyz = new double[grid.NumberOfBoundaryNodes][];
                double[][] boundaryUv = new double[grid.NumberOfBoundaryNodes][];
                int edge = 0;
                for (int l = 0; l < loopUv.Length; l++)
                {
                    for (int k = 0; k < loopUv[l].Length; k++, edge++)
                    {
                        if (poleTags.Contains(loopTags[l][k][0])) poleNodes.Add(grid.GetBoundaryNodeId(edge, 0));
                        for (int i = 0; i < grid.EdgeCounts[edge]; i++)
                        {
                            boundaryXyz[grid.GetBoundaryNodeId(edge, i)] = loopXyz[l][k][i];
                            boundaryUv[grid.GetBoundaryNodeId(edge, i)] = loopUv[l][k][i];
                        }
                    }
                }
                innerXyz = (double[][])innerXyz.Clone();
                innerUv = (double[][])innerUv.Clone();
                GmshQuadSurface.SmoothGrid(surfaceId, boundaryXyz, boundaryUv, innerXyz, innerUv, quads, poleNodes,
                                           ref numOfFoldedQuads, ref minScaledJacobian);
            }
            if (numOfFoldedQuads > 0)
            {
                reason = numOfFoldedQuads + " of the " + quads.Length + " quads of its grid fold";
                return false;
            }
            if (minScaledJacobian < SubmapFaceDetector.MinScaledJacobian)
            {
                reason = "its grid has a quad corner of the scaled Jacobian " + Math.Round(minScaledJacobian, 3) +
                         ", below " + SubmapFaceDetector.MinScaledJacobian;
                return false;
            }
            return true;
        }
        /// <summary>
        /// The number of elements a round hole of the radius aims at, given the number its edges are meshed with and
        /// the room beside it: with a ring around the hole the grid meets the circle a ring thickness outside it,
        /// which is longer, so the hole takes as many elements as keep those of that circle at elementSize, the
        /// length of the elements the grid beside the hole is made of. With sized = 2 pi radius / elementSize the
        /// ring, as thick as an element of the hole is long, gives n (n - sized) = 2 pi sized, and one of half the
        /// room, the thickest there is, n = sized (1 + room / (2 radius)). Never fewer than the hole is meshed with,
        /// which its curvature may have raised already, and that number where the room is less than MinRingRoom of
        /// elementSize: the rows more elements on the hole ask of the grid are then worse than the longer elements.
        /// </summary>
        public static double GetRingedNumberOfElements(int numOfElements, double radius, double room,
                                                       double elementSize)
        {
            if (!(elementSize > 0) || room < MinRingRoom * elementSize) return numOfElements;
            double sized = 2 * Math.PI * radius / elementSize;
            double ringed = 0.5 * (sized + Math.Sqrt(sized * sized + 8 * Math.PI * sized));
            if (2 * Math.PI * radius / ringed > 0.5 * room) ringed = sized * (1 + 0.5 * room / radius);
            return Math.Max(numOfElements, ringed);
        }
        /// <summary>
        /// The nodes the grid of the face is laid on: those of its loops, and for a round hole that is meshed with a
        /// ring those of the circle a ring thickness outside it, each on the line from the centre of the hole through
        /// its node of the hole and put on the surface.
        /// </summary>
        private static void GetRingLoops(int surfaceId, SubmapFace[] loops, double[][][][] loopXyz,
                                         double[][][][] loopUv, out double[][][][] gridXyz,
                                         out double[][][][] gridUv)
        {
            int numOfLoops = loops.Length;
            gridXyz = (double[][][][])loopXyz.Clone();
            gridUv = (double[][][][])loopUv.Clone();
            // The centre of every round hole, its radius and the length of its elements
            double[][] centres = new double[numOfLoops][];
            double[] radii = new double[numOfLoops];
            double[] lengths = new double[numOfLoops];
            for (int l = 1; l < numOfLoops; l++)
            {
                if (!loops[l].IsRound) continue;
                double[][] nodes = loopXyz[l].SelectMany(edge => edge.Take(edge.Length - 1)).ToArray();
                centres[l] = new double[3];
                foreach (var node in nodes)
                {
                    for (int d = 0; d < 3; d++) centres[l][d] += node[d] / nodes.Length;
                }
                for (int n = 0; n < nodes.Length; n++)
                {
                    radii[l] += GmshPoleChart.GetDistance(centres[l], nodes[n]) / nodes.Length;
                    lengths[l] += GmshPoleChart.GetDistance(nodes[n], nodes[(n + 1) % nodes.Length]) / nodes.Length;
                }
            }
            for (int l = 1; l < numOfLoops; l++)
            {
                if (centres[l] == null) continue;
                double room = double.MaxValue;
                for (int other = 0; other < numOfLoops; other++)
                {
                    if (other == l) continue;
                    double share = centres[other] != null ? 2.0 / 3 : 1;
                    foreach (var edge in loopXyz[other])
                    {
                        foreach (var node in edge)
                            room = Math.Min(room, share * (GmshPoleChart.GetDistance(centres[l], node) - radii[l]));
                    }
                }
                if (room < MinRingRoom * lengths[l]) continue;
                double thickness = Math.Min(lengths[l], 0.5 * room);
                int numOfEdges = loopXyz[l].Length;
                double[][][] ringXyz = new double[numOfEdges][][];
                double[][][] ringUv = new double[numOfEdges][][];
                bool placed = true;
                for (int k = 0; k < numOfEdges && placed; k++)
                {
                    int numOfNodes = loopXyz[l][k].Length;
                    ringXyz[k] = new double[numOfNodes][];
                    ringUv[k] = new double[numOfNodes][];
                    for (int i = 0; i < numOfNodes && placed; i++)
                    {
                        if (i == 0 && k > 0)
                        {
                            ringXyz[k][0] = ringXyz[k - 1][ringXyz[k - 1].Length - 1];
                            ringUv[k][0] = ringUv[k - 1][ringUv[k - 1].Length - 1];
                            continue;
                        }
                        if (i == numOfNodes - 1 && k == numOfEdges - 1)
                        {
                            ringXyz[k][i] = ringXyz[0][0];
                            ringUv[k][i] = ringUv[0][0];
                            continue;
                        }
                        double[] node = loopXyz[l][k][i];
                        double distance = GmshPoleChart.GetDistance(centres[l], node);
                        double[] point = new double[3];
                        for (int d = 0; d < 3; d++)
                            point[d] = node[d] + thickness * (node[d] - centres[l][d]) / distance;
                        double[] closest;
                        double[] parameters;
                        Gmsh.Model.GetClosestPoint(2, surfaceId, point, out closest, out _);
                        placed = closest.Length == 3;
                        if (!placed) break;
                        Gmsh.Model.GetParametrization(2, surfaceId, closest, out parameters);
                        placed = parameters.Length == 2;
                        if (!placed) break;
                        Gmsh.Model.GetValue(2, surfaceId, parameters, out closest);
                        ringXyz[k][i] = closest;
                        ringUv[k][i] = parameters;
                    }
                }
                if (!placed) continue;
                gridXyz[l] = ringXyz;
                gridUv[l] = ringUv;
            }
        }
        /// <summary>
        /// The quads and the inner nodes of the grid with the rings of its round holes: where the grid was laid on
        /// other nodes than those of a loop, gridUv against loopUv, its quads are joined to those nodes, which follow
        /// its inner nodes, and a ring of quads, turning as the quads of the grid, joins them to the nodes of the loop.
        /// The inner nodes, those the grid was laid on among them, are then moved together by the condition energy
        /// descent SubmapMapping ends with, in the frame of the outer loop, so that the five elements at a corner of
        /// the square of a hole share its angle. The number of folded quads and the worst scaled Jacobian of a quad
        /// corner are then those of all the quads in that frame; without a ring they are left as given.
        /// </summary>
        private static void AddRings(int surfaceId, SubmapGrid grid, double[][][][] loopXyz, double[][][][] loopUv,
                                     double[][][][] gridXyz, double[][][][] gridUv, out int[][] quads,
                                     out double[][] innerUv, out double[][] innerXyz, ref int numOfFoldedQuads,
                                     ref double minScaledJacobian)
        {
            quads = grid.Quads;
            innerUv = grid.InnerUv;
            innerXyz = grid.InnerXyz;
            if (Enumerable.Range(0, loopUv.Length).All(l => ReferenceEquals(loopUv[l], gridUv[l]))) return;
            List<double[]> allInnerUv = new List<double[]>(innerUv);
            List<double[]> allInnerXyz = new List<double[]>(innerXyz);
            Dictionary<int, int> ringNodeIds = new Dictionary<int, int>();
            List<int[]> ringSides = new List<int[]>();
            int edge = 0;
            for (int l = 0; l < loopUv.Length; l++)
            {
                for (int k = 0; k < loopUv[l].Length; k++, edge++)
                {
                    if (ReferenceEquals(loopUv[l], gridUv[l])) continue;
                    for (int i = 0; i < grid.EdgeCounts[edge]; i++)
                    {
                        ringNodeIds.Add(grid.GetBoundaryNodeId(edge, i), grid.NumberOfBoundaryNodes + allInnerUv.Count);
                        allInnerUv.Add(gridUv[l][k][i]);
                        allInnerXyz.Add(gridXyz[l][k][i]);
                        ringSides.Add(new int[] { grid.GetBoundaryNodeId(edge, i),
                                                  grid.GetBoundaryNodeId(edge, i + 1) });
                    }
                }
            }
            HashSet<Tuple<int, int>> gridSides = new HashSet<Tuple<int, int>>();
            List<int[]> allQuads = new List<int[]>();
            int ringNodeId;
            foreach (var quad in grid.Quads)
            {
                for (int c = 0; c < 4; c++) gridSides.Add(new Tuple<int, int>(quad[c], quad[(c + 1) % 4]));
                int[] joined = new int[4];
                for (int c = 0; c < 4; c++)
                    joined[c] = ringNodeIds.TryGetValue(quad[c], out ringNodeId) ? ringNodeId : quad[c];
                allQuads.Add(joined);
            }
            foreach (var side in ringSides)
            {
                int first = side[0];
                int second = side[1];
                if (!gridSides.Contains(new Tuple<int, int>(first, second)))
                {
                    first = side[1];
                    second = side[0];
                }
                allQuads.Add(new int[] { ringNodeIds[second], ringNodeIds[first], first, second });
            }
            quads = allQuads.ToArray();
            innerUv = allInnerUv.ToArray();
            // The frame of the outer loop, in which lengths and angles are those in space
            int numOfBoundaryNodes = grid.NumberOfBoundaryNodes;
            int numOfOuterNodes = grid.NumberOfOuterBoundaryNodes;
            int numOfNodes = numOfBoundaryNodes + innerUv.Length;
            double[][] boundaryUv = new double[numOfBoundaryNodes][];
            double[][] boundaryXyz = new double[numOfBoundaryNodes][];
            edge = 0;
            for (int l = 0; l < loopUv.Length; l++)
            {
                for (int k = 0; k < loopUv[l].Length; k++, edge++)
                {
                    for (int i = 0; i < grid.EdgeCounts[edge]; i++)
                    {
                        boundaryUv[grid.GetBoundaryNodeId(edge, i)] = loopUv[l][k][i];
                        boundaryXyz[grid.GetBoundaryNodeId(edge, i)] = loopXyz[l][k][i];
                    }
                }
            }
            double[][] outerUv = new double[numOfOuterNodes + 1][];
            double[][] outerXyz = new double[numOfOuterNodes + 1][];
            for (int b = 0; b <= numOfOuterNodes; b++)
            {
                outerUv[b] = boundaryUv[b % numOfOuterNodes];
                outerXyz[b] = boundaryXyz[b % numOfOuterNodes];
            }
            double[][] none = new double[0][];
            double[] w = ThreeBlockMapping.GetMetricFrame(new double[][][] { outerUv, none, none },
                                                          new double[][][] { outerXyz, none, none });
            double[][] points = new double[numOfNodes][];
            for (int b = 0; b < numOfBoundaryNodes; b++) points[b] = ThreeBlockMapping.ToFrame(w, boundaryUv[b]);
            for (int n = 0; n < innerUv.Length; n++)
                points[numOfBoundaryNodes + n] = ThreeBlockMapping.ToFrame(w, innerUv[n]);
            double orientation = ThreeBlockMapping.GetSignedArea(points, numOfOuterNodes) < 0 ? -1 : 1;
            int[][] neighbours = ThreeBlockMapping.GetNeighbours(numOfNodes, allQuads);
            int[] cornerStart;
            int[] cornerQuads;
            int[] cornerIndices;
            QuadSmoothing.GetCorners(quads, numOfNodes, numOfBoundaryNodes, out cornerStart, out cornerQuads,
                                     out cornerIndices);
            if (QuadSmoothing.CountFoldedElements(points, quads, orientation) > 0)
                QuadSmoothing.Untangle(points, neighbours, quads, cornerStart, cornerQuads, cornerIndices,
                                       numOfBoundaryNodes, orientation);
            QuadSmoothing.LowerConditionEnergy(points, neighbours, quads, cornerStart, cornerQuads, cornerIndices,
                                               numOfBoundaryNodes, orientation);
            numOfFoldedQuads = QuadSmoothing.CountFoldedElements(points, quads, orientation);
            minScaledJacobian = QuadSmoothing.GetMinScaledJacobian(points, quads, orientation);
            double[] flatUv = new double[2 * innerUv.Length];
            for (int n = 0; n < innerUv.Length; n++)
            {
                innerUv[n] = ThreeBlockMapping.FromFrame(w, points[numOfBoundaryNodes + n]);
                flatUv[2 * n] = innerUv[n][0];
                flatUv[2 * n + 1] = innerUv[n][1];
            }
            double[] flatXyz;
            Gmsh.Model.GetValue(2, surfaceId, flatUv, out flatXyz);
            innerXyz = new double[innerUv.Length][];
            for (int n = 0; n < innerUv.Length; n++)
                innerXyz[n] = new double[] { flatXyz[3 * n], flatXyz[3 * n + 1], flatXyz[3 * n + 2] };
        }
        internal static bool TryPlaceHoles(SubmapFace face, out List<IntPtr>[][] loopTags, out double[][][][] loopXyz,
                                           out double[][][][] loopUv, out double[][][] loopNodeUv,
                                           out int[][] holeOffsets, out string reason)
        {
            reason = null;
            holeOffsets = null;
            int surfaceId = face.FaceId;
            SubmapFace[] loops = face.Loops;
            int numOfLoops = loops.Length;
            int[][] loopDirections = face.LoopLogicalEdgeDirections;
            loopTags = new List<IntPtr>[numOfLoops][];
            loopXyz = new double[numOfLoops][][][];
            loopUv = new double[numOfLoops][][][];
            loopNodeUv = new double[numOfLoops][][];
            int[][] loopCounts = new int[numOfLoops][];
            for (int l = 0; l < numOfLoops; l++)
            {
                List<List<IntPtr>> logicalTags;
                List<double[][]> logicalXyz;
                List<double[][]> logicalUv;
                if (!GetLogicalEdgeNodes(surfaceId, loops[l], out logicalTags, out logicalXyz, out logicalUv, out reason))
                    return false;
                loopTags[l] = logicalTags.ToArray();
                loopXyz[l] = logicalXyz.ToArray();
                loopUv[l] = logicalUv.ToArray();
                loopCounts[l] = new int[logicalUv.Count];
                List<double[]> nodeUv = new List<double[]>();
                for (int k = 0; k < logicalUv.Count; k++)
                {
                    loopCounts[l][k] = logicalUv[k].Length - 1;
                    for (int i = 0; i < loopCounts[l][k]; i++) nodeUv.Add(logicalUv[k][i]);
                }
                loopNodeUv[l] = nodeUv.ToArray();
            }
            // A round hole starts its logical loop at the node the lattice of the outer loop puts nearest to its +i +j
            // corner
            for (int l = 1; l < numOfLoops; l++)
            {
                if (!loops[l].IsRound) continue;
                List<IntPtr>[] roundTags;
                double[][][] roundXyz;
                double[][][] roundUv;
                if (!GetRoundLoopEdgeNodes(loops[l], loopDirections[0], loopCounts[0], loopNodeUv[0], loopTags[l],
                                           loopXyz[l], loopUv[l], out roundTags, out roundXyz, out roundUv,
                                           out loopNodeUv[l], out reason)) return false;
                loopTags[l] = roundTags;
                loopXyz[l] = roundXyz;
                loopUv[l] = roundUv;
                loopCounts[l] = new int[4];
                for (int k = 0; k < 4; k++) loopCounts[l][k] = loopUv[l][k].Length - 1;
            }
            // The holes are placed in the lattice by where they lie in the parameter space
            if (!SubmapHolePlacement.TryPlace(loopDirections, loopCounts, loopNodeUv, 1, out holeOffsets, out reason))
            {
                reason = NameLoopEdges(face, reason);
                return false;
            }
            return true;
        }
        /// <summary>
        /// The nodes gmsh meshed on the logical edges of one loop of a face, by their tags, coordinates and parameters
        /// on the surface. A turn edge is split at its middle node into its two halves. False, with the reason, when the
        /// nodes cannot be read or a turn edge carries an odd number of elements.
        /// </summary>
        private static bool GetLogicalEdgeNodes(int surfaceId, SubmapFace loop, out List<List<IntPtr>> logicalTags,
                                                out List<double[][]> logicalXyz, out List<double[][]> logicalUv,
                                                out string reason)
        {
            reason = null;
            logicalTags = new List<List<IntPtr>>();
            logicalXyz = new List<double[][]>();
            logicalUv = new List<double[][]>();
            int numOfEdges = loop.EdgeIds.Length;
            List<IntPtr>[] edgeTags = new List<IntPtr>[numOfEdges];
            double[][][] edgeXyz = new double[numOfEdges][][];
            double[][][] edgeUv = new double[numOfEdges][][];
            try
            {
                for (int k = 0; k < numOfEdges; k++)
                {
                    List<double[]> xyz;
                    GmshThreeBlockMesher.GetSideNodes(surfaceId, loop.VertexNodeIds[k],
                                                      loop.VertexNodeIds[(k + 1) % numOfEdges],
                                                      new int[] { loop.EdgeIds[k] }, out edgeTags[k], out xyz);
                    edgeXyz[k] = xyz.ToArray();
                    edgeUv[k] = GmshThreeBlockMesher.GetParametricCoordinates(surfaceId, xyz);
                }
            }
            catch (CaeException ex)
            {
                reason = ToReason(ex.Message);
                return false;
            }
            // A turn edge is split at its middle node into its two halves, and an edge that turns twice at the nodes at a
            // quarter and three quarters of it into three pieces
            for (int k = 0; k < numOfEdges; k++)
            {
                int turn = Math.Abs(loop.GetEdgeTurn(k));
                if (turn == 0)
                {
                    logicalTags.Add(edgeTags[k]);
                    logicalXyz.Add(edgeXyz[k]);
                    logicalUv.Add(edgeUv[k]);
                    continue;
                }
                int count = edgeXyz[k].Length - 1;
                if (turn == 1 && count % 2 != 0)
                {
                    reason = "Edge " + loop.EdgeIds[k] + " turns a corner in its middle and carries an odd number of " +
                             "elements";
                    return false;
                }
                if (turn == 2 && count % 4 != 0)
                {
                    reason = "Edge " + loop.EdgeIds[k] + " turns twice and carries " + count + " elements, not a " +
                             "multiple of four";
                    return false;
                }
                int[] splits = turn == 1 ? new int[] { 0, count / 2, count } :
                                           new int[] { 0, count / 4, 3 * count / 4, count };
                for (int p = 0; p < splits.Length - 1; p++)
                {
                    int length = splits[p + 1] - splits[p] + 1;
                    logicalTags.Add(edgeTags[k].GetRange(splits[p], length));
                    logicalXyz.Add(GetRange(edgeXyz[k], splits[p], length));
                    logicalUv.Add(GetRange(edgeUv[k], splits[p], length));
                }
            }
            return true;
        }
        /// <summary>
        /// The four logical edges of a round hole from the nodes gmsh meshed on its edges, which must carry as many
        /// elements each and a multiple of four around it: its nodes in the order SubmapHolePlacement.GetRoundLoopOrder
        /// finds in the lattice of the outer loop, whose logical edges run in outerDirections with outerCounts elements
        /// through the nodes of the parameters outerNodeUv, and nodeUv the parameters of its nodes in that order. False,
        /// with the reason, when the counts do not allow it.
        /// </summary>
        private static bool GetRoundLoopEdgeNodes(SubmapFace loop, int[] outerDirections, int[] outerCounts,
                                                  double[][] outerNodeUv, List<IntPtr>[] edgeTags, double[][][] edgeXyz,
                                                  double[][][] edgeUv, out List<IntPtr>[] logicalTags,
                                                  out double[][][] logicalXyz, out double[][][] logicalUv,
                                                  out double[][] nodeUv, out string reason)
        {
            logicalTags = null;
            logicalXyz = null;
            logicalUv = null;
            nodeUv = null;
            reason = null;
            List<IntPtr> ringTags = new List<IntPtr>();
            List<double[]> ringXyz = new List<double[]>();
            List<double[]> ringUv = new List<double[]>();
            int firstCount = edgeUv[0].Length - 1;
            for (int k = 0; k < edgeUv.Length; k++)
            {
                int count = edgeUv[k].Length - 1;
                if (count != firstCount && !loop.HasUnequalRoundEdges)
                {
                    reason = "Edge " + loop.EdgeIds[0] + " and Edge " + loop.EdgeIds[k] + " of a round hole carry " +
                             firstCount + " and " + count + " elements, which must be equal";
                    return false;
                }
                for (int i = 0; i < count; i++)
                {
                    ringTags.Add(edgeTags[k][i]);
                    ringXyz.Add(edgeXyz[k][i]);
                    ringUv.Add(edgeUv[k][i]);
                }
            }
            int numOfNodes = ringUv.Count;
            if (numOfNodes % 4 != 0 || numOfNodes == 0)
            {
                reason = "the round hole of Edge " + loop.EdgeIds[0] + " carries " + numOfNodes + " elements around it, " +
                         "not a multiple of four";
                return false;
            }
            int[] order = SubmapHolePlacement.GetRoundLoopOrder(outerDirections, outerCounts, outerNodeUv,
                                                                ringUv.ToArray());
            int edgeCount = numOfNodes / 4;
            nodeUv = new double[numOfNodes][];
            for (int p = 0; p < numOfNodes; p++) nodeUv[p] = ringUv[order[p]];
            logicalTags = new List<IntPtr>[4];
            logicalXyz = new double[4][][];
            logicalUv = new double[4][][];
            for (int c = 0; c < 4; c++)
            {
                logicalTags[c] = new List<IntPtr>();
                logicalXyz[c] = new double[edgeCount + 1][];
                logicalUv[c] = new double[edgeCount + 1][];
                for (int i = 0; i <= edgeCount; i++)
                {
                    int node = order[(c * edgeCount + i) % numOfNodes];
                    logicalTags[c].Add(ringTags[node]);
                    logicalXyz[c][i] = ringXyz[node];
                    logicalUv[c][i] = ringUv[node];
                }
            }
            return true;
        }
        /// <summary>
        /// A reason SubmapMapping gives, which names the edges of the face by their places in all its loops one after
        /// another, with the edges named by their gmsh ids, as a clause that follows a dash.
        /// </summary>
        internal static string NameLoopEdges(SubmapFace face, string reason)
        {
            reason = LoopEdge.Replace(reason, match =>
            {
                string text = match.Groups[1].Value.Length > 0 ? "Edges " : "Edge ";
                text += GetEdgeId(face, match.Groups[2].Value);
                if (match.Groups[3].Success) text += " and " + GetEdgeId(face, match.Groups[3].Value);
                return text;
            });
            return ToReason(reason);
        }
        private static string GetEdgeId(SubmapFace face, string loopIndex)
        {
            int index = int.Parse(loopIndex);
            int[] logicalEdgeIds = face.AllLogicalEdgeIds;
            return index >= 0 && index < logicalEdgeIds.Length ? logicalEdgeIds[index].ToString() : loopIndex;
        }
        private static double[][] GetRange(double[][] array, int start, int length)
        {
            double[][] range = new double[length][];
            Array.Copy(array, start, range, 0, length);
            return range;
        }
        private static string ToReason(string message)
        {
            message = message.Trim().TrimEnd('.');
            if (message.StartsWith("The ")) message = "the " + message.Substring(4);
            return message;
        }
    }
}
