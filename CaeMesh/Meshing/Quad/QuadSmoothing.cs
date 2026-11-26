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

namespace CaeMesh
{
    /// <summary>
    /// The smoothing every quad mesher shares, in a plane frame close to isometric, the boundary nodes fixed. Its measure
    /// of a corner is the scaled Jacobian, the sine of the corner angle, negative where the element turns against the
    /// orientation of the face, and for a triangle the sine over that of 60 degrees, so that an equilateral triangle is
    /// as good as a square. The smoothing lowers the condition energy of the corners a node changes: (|a|^2 + |b|^2) /
    /// (2 a x b) at a quad corner and (|a|^2 + |b|^2 - a.b) / (sqrt(3) a x b) at a triangle corner, a and b the edges at
    /// the corner, 1 at a square or an equilateral corner, growing with skew and stretch and infinite for a corner that
    /// folds, so that no step folds an element. A node is moved by a pattern search over eight directions with a step
    /// halved six times. Folded elements are first untangled by moves that raise the worst corner of a node.
    /// The structured grids number their boundary nodes first and smooth every node after them (Untangle,
    /// LowerConditionEnergy, RelaxWithoutFolding). An unstructured region of quads and triangles of any valence is
    /// smoothed on a given set of free nodes by SmoothRegion, which also never lowers the worst corner of the elements
    /// it touches. SmoothRegionInSpace smooths such a region on a surface in space, for which no plane frame is needed.
    /// The arithmetic is plain and ordered, so that a port gives the same mesh.
    /// </summary>
    public static class QuadSmoothing
    {
        // Variables                                                                                                                
        private const int MaxNumOfUntangleSweeps = 100;
        private const int NumOfUntangleStepTrials = 8;
        private const int MaxNumOfDescentSweeps = 200;
        private const int NumOfDescentStepTrials = 6;
        private const double DescentStepFactor = 0.1;
        private const double DescentMinImprovement = 1E-9;
        private const int MaxNumOfGuardedSweeps = 500;
        private const int NumOfGuardedStepTrials = 4;
        private const double GuardedTolerance = 1E-10;
        private const double RegionMoveTolerance = 1E-3;
        private const double RegionMinRaise = 1E-4;
        private const double SlopeStepFactor = 1E-3;
        private static readonly int[] CornerOffsets = new int[] { 0, 1, 3 };
        private static readonly double TriangleScale = 2 / Math.Sqrt(3);
        private static readonly double SqrtThree = Math.Sqrt(3);


        // Methods                                                                                                                  
        /// <summary>
        /// The corners every inner node changes, as quad index and corner index, the nodes after each other: of every
        /// quad, every corner in turn, and of a corner its node and the nodes after it and before it.
        /// </summary>
        internal static void GetCorners(int[][] quads, int numOfNodes, int numOfBoundaryNodes, out int[] cornerStart,
                                       out int[] cornerQuads, out int[] cornerIndices)
        {
            cornerStart = new int[numOfNodes + 1];
            foreach (var quad in quads)
            {
                for (int c = 0; c < 4; c++)
                {
                    for (int o = 0; o < 3; o++)
                    {
                        int node = quad[(c + CornerOffsets[o]) % 4];
                        if (node >= numOfBoundaryNodes) cornerStart[node + 1]++;
                    }
                }
            }
            for (int n = 0; n < numOfNodes; n++) cornerStart[n + 1] += cornerStart[n];
            cornerQuads = new int[cornerStart[numOfNodes]];
            cornerIndices = new int[cornerQuads.Length];
            int[] fill = (int[])cornerStart.Clone();
            for (int q = 0; q < quads.Length; q++)
            {
                for (int c = 0; c < 4; c++)
                {
                    for (int o = 0; o < 3; o++)
                    {
                        int node = quads[q][(c + CornerOffsets[o]) % 4];
                        if (node < numOfBoundaryNodes) continue;
                        cornerQuads[fill[node]] = q;
                        cornerIndices[fill[node]] = c;
                        fill[node]++;
                    }
                }
            }
        }
        /// <summary>
        /// Moves the inner nodes of folded quads, in turn and sweep after sweep, by a pattern search that only takes a
        /// step when it raises the smallest scaled Jacobian of the quad corners the node changes, until no quad folds
        /// or no node moves.
        /// </summary>
        internal static void Untangle(double[][] points, int[][] neighbours, int[][] quads, int[] cornerStart,
                                      int[] cornerQuads, int[] cornerIndices, int numOfBoundaryNodes, double orientation)
        {
            int numOfNodes = points.Length;
            double[][] directions = GetSearchDirections();
            //
            for (int sweep = 0; sweep < MaxNumOfUntangleSweeps; sweep++)
            {
                bool folded = false;
                bool moved = false;
                for (int node = numOfBoundaryNodes; node < numOfNodes; node++)
                {
                    double current = GetMinScaledJacobian(points, quads, cornerStart, cornerQuads, cornerIndices, node,
                                                          orientation);
                    if (current > 0) continue;
                    folded = true;
                    double[] point = points[node];
                    double length = 0;
                    foreach (var neighbour in neighbours[node])
                    {
                        double dx = points[neighbour][0] - point[0];
                        double dy = points[neighbour][1] - point[1];
                        length += Math.Sqrt(dx * dx + dy * dy);
                    }
                    double step = 0.5 * length / neighbours[node].Length;
                    for (int trial = 0; trial < NumOfUntangleStepTrials; trial++)
                    {
                        double x0 = point[0];
                        double y0 = point[1];
                        double best = current + 1E-12;
                        int bestDirection = -1;
                        for (int d = 0; d < directions.Length; d++)
                        {
                            point[0] = x0 + directions[d][0] * step;
                            point[1] = y0 + directions[d][1] * step;
                            double q = GetMinScaledJacobian(points, quads, cornerStart, cornerQuads, cornerIndices, node,
                                                            orientation);
                            if (q > best)
                            {
                                best = q;
                                bestDirection = d;
                            }
                        }
                        if (bestDirection >= 0)
                        {
                            point[0] = x0 + directions[bestDirection][0] * step;
                            point[1] = y0 + directions[bestDirection][1] * step;
                            current = best;
                            moved = true;
                        }
                        else
                        {
                            point[0] = x0;
                            point[1] = y0;
                            step *= 0.5;
                        }
                    }
                }
                if (!folded || !moved) break;
            }
        }
        /// <summary>
        /// The condition energy descent over the inner nodes of a quad grid, run until no node moves.
        /// </summary>
        internal static void LowerConditionEnergy(double[][] points, int[][] neighbours, int[][] quads, int[] cornerStart,
                                                  int[] cornerQuads, int[] cornerIndices, int numOfBoundaryNodes,
                                                  double orientation)
        {
            LowerConditionEnergy(points, neighbours, quads, cornerStart, cornerQuads, cornerIndices, numOfBoundaryNodes,
                                 orientation, MaxNumOfDescentSweeps);
        }
        /// <summary>
        /// Moves the inner nodes, in turn and sweep after sweep, by a pattern search that only takes a step when it lowers
        /// the condition energy of the quad corners the node changes, for at most the given number of sweeps. A node is
        /// visited only when it or a node of its quads moved since its last visit; any other visit would find the same
        /// energies and move nothing, so the sweeps give the same grid as sweeps over every node, but a sweep over a
        /// large grid that has settled costs little.
        /// </summary>
        internal static void LowerConditionEnergy(double[][] points, int[][] neighbours, int[][] quads, int[] cornerStart,
                                                  int[] cornerQuads, int[] cornerIndices, int numOfBoundaryNodes,
                                                  double orientation, int maxNumOfSweeps)
        {
            int numOfNodes = points.Length;
            bool[] changed = new bool[numOfNodes];
            for (int node = numOfBoundaryNodes; node < numOfNodes; node++) changed[node] = true;
            double[][] directions = GetSearchDirections();
            for (int sweep = 0; sweep < maxNumOfSweeps; sweep++)
            {
                bool moved = false;
                for (int node = numOfBoundaryNodes; node < numOfNodes; node++)
                {
                    if (!changed[node]) continue;
                    changed[node] = false;
                    if (cornerStart[node] == cornerStart[node + 1]) continue;
                    double[] point = points[node];
                    double length = 0;
                    foreach (var neighbour in neighbours[node])
                    {
                        double dx = points[neighbour][0] - point[0];
                        double dy = points[neighbour][1] - point[1];
                        length += Math.Sqrt(dx * dx + dy * dy);
                    }
                    double step = DescentStepFactor * length / neighbours[node].Length;
                    double current = GetNodeQuality(points, quads, cornerStart, cornerQuads, cornerIndices, node,
                                                    orientation);
                    bool nodeMoved = false;
                    for (int trial = 0; trial < NumOfDescentStepTrials; trial++)
                    {
                        double x0 = point[0];
                        double y0 = point[1];
                        double bestQuality = double.IsNegativeInfinity(current) ? double.MinValue :
                                             current + DescentMinImprovement * Math.Max(1, Math.Abs(current));
                        int best = -1;
                        for (int d = 0; d < directions.Length; d++)
                        {
                            point[0] = x0 + directions[d][0] * step;
                            point[1] = y0 + directions[d][1] * step;
                            double q = GetNodeQuality(points, quads, cornerStart, cornerQuads, cornerIndices, node,
                                                      orientation);
                            if (q > bestQuality)
                            {
                                bestQuality = q;
                                best = d;
                            }
                        }
                        if (best >= 0)
                        {
                            point[0] = x0 + directions[best][0] * step;
                            point[1] = y0 + directions[best][1] * step;
                            current = bestQuality;
                            nodeMoved = true;
                        }
                        else
                        {
                            point[0] = x0;
                            point[1] = y0;
                            step *= 0.5;
                        }
                    }
                    if (!nodeMoved) continue;
                    moved = true;
                    for (int i = cornerStart[node]; i < cornerStart[node + 1]; i++)
                    {
                        foreach (var quadNode in quads[cornerQuads[i]]) changed[quadNode] = true;
                    }
                }
                if (!moved) break;
            }
        }
        /// <summary>
        /// Moves the nodes after the boundary ones, in turn and sweep after sweep, towards the mean of their neighbours,
        /// by the whole way or by a half, a quarter or an eighth of it, whichever first folds none of their quads, until
        /// no node moves by more than the tolerance; a grid that does not fold stays so.
        /// </summary>
        internal static void RelaxWithoutFolding(double[][] points, int[][] neighbours, int[][] quads,
                                                 int numOfBoundaryNodes, double orientation)
        {
            int numOfNodes = points.Length;
            int[] quadStart = new int[numOfNodes + 1];
            foreach (var quad in quads)
            {
                for (int c = 0; c < 4; c++) quadStart[quad[c] + 1]++;
            }
            for (int n = 0; n < numOfNodes; n++) quadStart[n + 1] += quadStart[n];
            int[] nodeQuads = new int[quadStart[numOfNodes]];
            int[] fill = (int[])quadStart.Clone();
            for (int q = 0; q < quads.Length; q++)
            {
                for (int c = 0; c < 4; c++) nodeQuads[fill[quads[q][c]]++] = q;
            }
            double minX = double.MaxValue;
            double maxX = double.MinValue;
            double minY = double.MaxValue;
            double maxY = double.MinValue;
            for (int b = 0; b < numOfBoundaryNodes; b++)
            {
                minX = Math.Min(minX, points[b][0]);
                maxX = Math.Max(maxX, points[b][0]);
                minY = Math.Min(minY, points[b][1]);
                maxY = Math.Max(maxY, points[b][1]);
            }
            double tolerance = GuardedTolerance * Math.Max(maxX - minX, maxY - minY);
            //
            for (int sweep = 0; sweep < MaxNumOfGuardedSweeps; sweep++)
            {
                double maxMove = 0;
                for (int node = numOfBoundaryNodes; node < numOfNodes; node++)
                {
                    double[] point = points[node];
                    double x0 = point[0];
                    double y0 = point[1];
                    double meanX = 0;
                    double meanY = 0;
                    foreach (var neighbour in neighbours[node])
                    {
                        meanX += points[neighbour][0];
                        meanY += points[neighbour][1];
                    }
                    meanX /= neighbours[node].Length;
                    meanY /= neighbours[node].Length;
                    double fraction = 1;
                    for (int trial = 0; trial < NumOfGuardedStepTrials; trial++)
                    {
                        point[0] = x0 + (meanX - x0) * fraction;
                        point[1] = y0 + (meanY - y0) * fraction;
                        if (!IsAnyQuadFolded(points, quads, nodeQuads, quadStart[node], quadStart[node + 1], orientation))
                            break;
                        point[0] = x0;
                        point[1] = y0;
                        fraction *= 0.5;
                    }
                    maxMove = Math.Max(maxMove, Math.Max(Math.Abs(point[0] - x0), Math.Abs(point[1] - y0)));
                }
                if (maxMove <= tolerance) break;
            }
        }
        /// <summary>
        /// Smooths the given free nodes of a region of quads and triangles of any valence, every node of the elements
        /// given that is not in the list staying where it is. The elements must hold every element of a free node, and
        /// turn the way the orientation says. The nodes of a corner that folds or is below the target are first moved by
        /// a pattern search that raises the worst corner they change, and the condition energy is then lowered for at
        /// most the given number of sweeps, a step taken only when no corner it changes falls below the worst corner of
        /// the elements after the first stage, so that the worst corner never falls and no element that does not fold
        /// starts to. The mean of the neighbours of a node is tried first, and the pattern search only when that is no
        /// better; a node is visited again only when a node of its elements moved by more than a thousandth of its mean
        /// edge length. Returns the worst scaled Jacobian of the elements.
        /// </summary>
        public static double SmoothRegion(double[][] points, int[][] elements, int[] freeNodes, double orientation,
                                          int maxNumOfSweeps, double target)
        {
            int numOfFree = freeNodes.Length;
            Dictionary<int, int> nodeIndex;
            int[] cornerStart;
            int[] cornerElements;
            int[] cornerIndices;
            int[][] neighbours;
            GetRegionCorners(elements, freeNodes, out nodeIndex, out cornerStart, out cornerElements, out cornerIndices,
                             out neighbours);
            int index;
            // Untangling, and the worst corners raised
            double floor = GetMinScaledJacobian(points, elements, orientation);
            double[][] directions = GetSearchDirections();
            if (floor <= 0 || floor < target)
            {
                int maxNumOfRaiseSweeps = Math.Min(MaxNumOfUntangleSweeps, 3 * maxNumOfSweeps);
                for (int sweep = 0; sweep < maxNumOfRaiseSweeps; sweep++)
                {
                    bool folded = false;
                    bool moved = false;
                    double maxRaise = 0;
                    for (int i = 0; i < numOfFree; i++)
                    {
                        double energy;
                        double current = GetRegionNodeQuality(points, elements, cornerElements, cornerIndices,
                                                              cornerStart[i], cornerStart[i + 1], orientation, out energy);
                        if ((current > 0 && current >= target) || neighbours[i].Length == 0) continue;
                        folded = true;
                        double[] point = points[freeNodes[i]];
                        double step = 0.5 * GetMeanLength(points, point, neighbours[i]);
                        double start = current;
                        for (int trial = 0; trial < NumOfUntangleStepTrials; trial++)
                        {
                            double x0 = point[0];
                            double y0 = point[1];
                            double best = current + 1E-12;
                            int bestDirection = -1;
                            for (int d = 0; d < directions.Length; d++)
                            {
                                point[0] = x0 + directions[d][0] * step;
                                point[1] = y0 + directions[d][1] * step;
                                double q = GetRegionNodeQuality(points, elements, cornerElements, cornerIndices,
                                                                cornerStart[i], cornerStart[i + 1], orientation,
                                                                out energy);
                                if (q > best)
                                {
                                    best = q;
                                    bestDirection = d;
                                }
                            }
                            if (bestDirection >= 0)
                            {
                                point[0] = x0 + directions[bestDirection][0] * step;
                                point[1] = y0 + directions[bestDirection][1] * step;
                                current = best;
                                moved = true;
                            }
                            else
                            {
                                point[0] = x0;
                                point[1] = y0;
                                step *= 0.5;
                            }
                        }
                        maxRaise = Math.Max(maxRaise, current - start);
                    }
                    if (!folded || !moved || maxRaise < RegionMinRaise) break;
                }
                floor = GetMinScaledJacobian(points, elements, orientation);
            }
            // The condition energy descent, guarded by the worst corner
            bool[] changed = new bool[numOfFree];
            for (int i = 0; i < numOfFree; i++) changed[i] = true;
            for (int sweep = 0; sweep < maxNumOfSweeps; sweep++)
            {
                bool moved = false;
                for (int i = 0; i < numOfFree; i++)
                {
                    if (!changed[i]) continue;
                    changed[i] = false;
                    if (cornerStart[i] == cornerStart[i + 1] || neighbours[i].Length == 0) continue;
                    double[] point = points[freeNodes[i]];
                    double energy;
                    GetRegionNodeQuality(points, elements, cornerElements, cornerIndices, cornerStart[i],
                                         cornerStart[i + 1], orientation, out energy);
                    double current = -energy;
                    bool nodeMoved = false;
                    // The mean of the neighbours first
                    double startX = point[0];
                    double startY = point[1];
                    double x0 = point[0];
                    double y0 = point[1];
                    double meanX = 0;
                    double meanY = 0;
                    foreach (var neighbour in neighbours[i])
                    {
                        meanX += points[neighbour][0];
                        meanY += points[neighbour][1];
                    }
                    point[0] = meanX / neighbours[i].Length;
                    point[1] = meanY / neighbours[i].Length;
                    double threshold = double.IsNegativeInfinity(current) ? double.MinValue :
                                       current + DescentMinImprovement * Math.Max(1, Math.Abs(current));
                    double minCorner = GetRegionNodeQuality(points, elements, cornerElements, cornerIndices,
                                                            cornerStart[i], cornerStart[i + 1], orientation, out energy);
                    if (-energy > threshold && minCorner >= floor)
                    {
                        current = -energy;
                        nodeMoved = true;
                    }
                    else
                    {
                        point[0] = x0;
                        point[1] = y0;
                    }
                    double meanLength = GetMeanLength(points, point, neighbours[i]);
                    double step = DescentStepFactor * meanLength;
                    for (int trial = 0; trial < (nodeMoved ? 0 : NumOfDescentStepTrials); trial++)
                    {
                        x0 = point[0];
                        y0 = point[1];
                        double bestQuality = double.IsNegativeInfinity(current) ? double.MinValue :
                                             current + DescentMinImprovement * Math.Max(1, Math.Abs(current));
                        int best = -1;
                        for (int d = 0; d < directions.Length; d++)
                        {
                            point[0] = x0 + directions[d][0] * step;
                            point[1] = y0 + directions[d][1] * step;
                            minCorner = GetRegionNodeQuality(points, elements, cornerElements, cornerIndices,
                                                             cornerStart[i], cornerStart[i + 1], orientation, out energy);
                            if (-energy > bestQuality && minCorner >= floor)
                            {
                                bestQuality = -energy;
                                best = d;
                            }
                        }
                        if (best >= 0)
                        {
                            point[0] = x0 + directions[best][0] * step;
                            point[1] = y0 + directions[best][1] * step;
                            current = bestQuality;
                            nodeMoved = true;
                        }
                        else
                        {
                            point[0] = x0;
                            point[1] = y0;
                            step *= 0.5;
                        }
                    }
                    if (!nodeMoved) continue;
                    double dx = point[0] - startX;
                    double dy = point[1] - startY;
                    if (!(Math.Sqrt(dx * dx + dy * dy) > RegionMoveTolerance * meanLength)) continue;
                    moved = true;
                    for (int k = cornerStart[i]; k < cornerStart[i + 1]; k++)
                    {
                        foreach (var elementNode in elements[cornerElements[k]])
                        {
                            if (nodeIndex.TryGetValue(elementNode, out index)) changed[index] = true;
                        }
                    }
                }
                if (!moved) break;
            }
            return GetMinScaledJacobian(points, elements, orientation);
        }
        /// <summary>
        /// The corners every free node of a region changes, as element index and corner index, and its neighbours along
        /// the element edges.
        /// </summary>
        private static void GetRegionCorners(int[][] elements, int[] freeNodes, out Dictionary<int, int> nodeIndex,
                                             out int[] cornerStart, out int[] cornerElements, out int[] cornerIndices,
                                             out int[][] neighbours)
        {
            int numOfFree = freeNodes.Length;
            nodeIndex = new Dictionary<int, int>();
            for (int i = 0; i < numOfFree; i++) nodeIndex[freeNodes[i]] = i;
            // The corners every free node changes, and its neighbours along the element edges
            cornerStart = new int[numOfFree + 1];
            int index;
            foreach (var element in elements)
            {
                int n = element.Length;
                for (int p = 0; p < n; p++)
                {
                    if (nodeIndex.TryGetValue(element[p], out index)) cornerStart[index + 1] += 3;
                }
            }
            for (int i = 0; i < numOfFree; i++) cornerStart[i + 1] += cornerStart[i];
            cornerElements = new int[cornerStart[numOfFree]];
            cornerIndices = new int[cornerElements.Length];
            int[] fill = (int[])cornerStart.Clone();
            List<int>[] neighbourLists = new List<int>[numOfFree];
            for (int i = 0; i < numOfFree; i++) neighbourLists[i] = new List<int>();
            for (int e = 0; e < elements.Length; e++)
            {
                int[] element = elements[e];
                int n = element.Length;
                for (int p = 0; p < n; p++)
                {
                    if (!nodeIndex.TryGetValue(element[p], out index)) continue;
                    // A quad corner is changed by its node and the two next to it, a triangle corner by all three
                    int[] corners = new int[] { p, (p + n - 1) % n, (p + 1) % n };
                    foreach (var c in corners)
                    {
                        cornerElements[fill[index]] = e;
                        cornerIndices[fill[index]] = c;
                        fill[index]++;
                    }
                    AddUnique(neighbourLists[index], element[(p + 1) % n]);
                    AddUnique(neighbourLists[index], element[(p + n - 1) % n]);
                }
            }
            neighbours = new int[numOfFree][];
            for (int i = 0; i < numOfFree; i++)
            {
                neighbourLists[i].Sort();
                neighbours[i] = neighbourLists[i].ToArray();
            }
        }
        /// <summary>
        /// SmoothRegion on a surface in space, which needs no plane the whole region can be laid in: a corner is measured
        /// in space against the normal of the surface at its node times the orientation - at a node without a normal, a
        /// pole of the surface, against the normal of its own element - and a node moves in the tangent plane of the
        /// surface at it, first the way its measure falls fastest (TurnTangentsDownhill). The nodes a sweep moved are
        /// put back onto the surface when the sweep ends (IQuadSurface.Project), which also renews their normals: a
        /// move in the tangent plane leaves the surface by a length of the second order only. Returns the worst scaled
        /// Jacobian of the elements.
        /// </summary>
        public static double SmoothRegionInSpace(double[][] xyz, double[][] normals, int[][] elements, int[] freeNodes,
                                                 double orientation, int maxNumOfSweeps, double target,
                                                 IQuadSurface surface)
        {
            int numOfFree = freeNodes.Length;
            Dictionary<int, int> nodeIndex;
            int[] cornerStart;
            int[] cornerElements;
            int[] cornerIndices;
            int[][] neighbours;
            GetRegionCorners(elements, freeNodes, out nodeIndex, out cornerStart, out cornerElements, out cornerIndices,
                             out neighbours);
            int index;
            double energy;
            double[] e1 = new double[3];
            double[] e2 = new double[3];
            double[] origin = new double[3];
            List<int> movedNodes = new List<int>();
            List<double[]> previousXyz = new List<double[]>();
            // Untangling, and the worst corners raised
            double floor = GetSpaceMinScaledJacobian(xyz, normals, elements, orientation);
            double[][] directions = GetSearchDirections();
            if (floor <= 0 || floor < target)
            {
                int maxNumOfRaiseSweeps = Math.Min(MaxNumOfUntangleSweeps, 3 * maxNumOfSweeps);
                for (int sweep = 0; sweep < maxNumOfRaiseSweeps; sweep++)
                {
                    bool folded = false;
                    double maxRaise = 0;
                    for (int i = 0; i < numOfFree; i++)
                    {
                        int node = freeNodes[i];
                        if (normals[node] == null || neighbours[i].Length == 0) continue;
                        double current = GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements, cornerIndices,
                                                                   cornerStart[i], cornerStart[i + 1], orientation,
                                                                   out energy);
                        if (current > 0 && current >= target) continue;
                        folded = true;
                        double[] previous = xyz[node];
                        double[] point = (double[])previous.Clone();
                        xyz[node] = point;
                        double step = 0.5 * GetSpaceMeanLength(xyz, point, neighbours[i]);
                        GetTangents(normals[node], e1, e2);
                        TurnTangentsDownhill(xyz, normals, elements, cornerElements, cornerIndices, cornerStart[i],
                                             cornerStart[i + 1], orientation, point, normals[node],
                                             SlopeStepFactor * step, false, e1, e2);
                        double start = current;
                        bool nodeMoved = false;
                        for (int trial = 0; trial < NumOfUntangleStepTrials; trial++)
                        {
                            Array.Copy(point, origin, 3);
                            double best = current + 1E-12;
                            int bestDirection = -1;
                            for (int d = 0; d < directions.Length; d++)
                            {
                                SetTangentPoint(point, origin, e1, e2, directions[d][0] * step, directions[d][1] * step);
                                double q = GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements,
                                                                     cornerIndices, cornerStart[i], cornerStart[i + 1],
                                                                     orientation, out energy);
                                if (q > best)
                                {
                                    best = q;
                                    bestDirection = d;
                                }
                            }
                            if (bestDirection >= 0)
                            {
                                SetTangentPoint(point, origin, e1, e2, directions[bestDirection][0] * step,
                                                directions[bestDirection][1] * step);
                                current = best;
                                nodeMoved = true;
                            }
                            else
                            {
                                Array.Copy(origin, point, 3);
                                step *= 0.5;
                            }
                        }
                        if (nodeMoved)
                        {
                            movedNodes.Add(node);
                            previousXyz.Add(previous);
                        }
                        else xyz[node] = previous;
                        maxRaise = Math.Max(maxRaise, current - start);
                    }
                    bool moved = movedNodes.Count > 0;
                    surface.Project(xyz, normals, movedNodes, previousXyz.ToArray());
                    movedNodes.Clear();
                    previousXyz.Clear();
                    if (!folded || !moved || maxRaise < RegionMinRaise) break;
                }
                floor = GetSpaceMinScaledJacobian(xyz, normals, elements, orientation);
            }
            // The condition energy descent, guarded by the worst corner
            bool[] changed = new bool[numOfFree];
            for (int i = 0; i < numOfFree; i++) changed[i] = true;
            for (int sweep = 0; sweep < maxNumOfSweeps; sweep++)
            {
                bool moved = false;
                for (int i = 0; i < numOfFree; i++)
                {
                    if (!changed[i]) continue;
                    changed[i] = false;
                    int node = freeNodes[i];
                    if (cornerStart[i] == cornerStart[i + 1] || neighbours[i].Length == 0 || normals[node] == null)
                        continue;
                    double[] normal = normals[node];
                    double[] previous = xyz[node];
                    double[] point = (double[])previous.Clone();
                    xyz[node] = point;
                    GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements, cornerIndices, cornerStart[i],
                                              cornerStart[i + 1], orientation, out energy);
                    double current = -energy;
                    bool nodeMoved = false;
                    // The mean of the neighbours first, in the tangent plane
                    double[] mean = new double[3];
                    foreach (var neighbour in neighbours[i])
                    {
                        for (int d = 0; d < 3; d++) mean[d] += xyz[neighbour][d];
                    }
                    double height = 0;
                    for (int d = 0; d < 3; d++)
                    {
                        mean[d] /= neighbours[i].Length;
                        height += (mean[d] - previous[d]) * normal[d];
                    }
                    for (int d = 0; d < 3; d++) point[d] = mean[d] - height * normal[d];
                    double threshold = double.IsNegativeInfinity(current) ? double.MinValue :
                                       current + DescentMinImprovement * Math.Max(1, Math.Abs(current));
                    double minCorner = GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements, cornerIndices,
                                                                 cornerStart[i], cornerStart[i + 1], orientation,
                                                                 out energy);
                    if (-energy > threshold && minCorner >= floor)
                    {
                        current = -energy;
                        nodeMoved = true;
                    }
                    else Array.Copy(previous, point, 3);
                    double meanLength = GetSpaceMeanLength(xyz, point, neighbours[i]);
                    double step = DescentStepFactor * meanLength;
                    if (!nodeMoved)
                    {
                        GetTangents(normal, e1, e2);
                        TurnTangentsDownhill(xyz, normals, elements, cornerElements, cornerIndices, cornerStart[i],
                                             cornerStart[i + 1], orientation, point, normal, SlopeStepFactor * step,
                                             true, e1, e2);
                    }
                    for (int trial = 0; trial < (nodeMoved ? 0 : NumOfDescentStepTrials); trial++)
                    {
                        Array.Copy(point, origin, 3);
                        double bestQuality = double.IsNegativeInfinity(current) ? double.MinValue :
                                             current + DescentMinImprovement * Math.Max(1, Math.Abs(current));
                        int best = -1;
                        for (int d = 0; d < directions.Length; d++)
                        {
                            SetTangentPoint(point, origin, e1, e2, directions[d][0] * step, directions[d][1] * step);
                            minCorner = GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements, cornerIndices,
                                                                  cornerStart[i], cornerStart[i + 1], orientation,
                                                                  out energy);
                            if (-energy > bestQuality && minCorner >= floor)
                            {
                                bestQuality = -energy;
                                best = d;
                            }
                        }
                        if (best >= 0)
                        {
                            SetTangentPoint(point, origin, e1, e2, directions[best][0] * step, directions[best][1] * step);
                            current = bestQuality;
                            nodeMoved = true;
                        }
                        else
                        {
                            Array.Copy(origin, point, 3);
                            step *= 0.5;
                        }
                    }
                    if (!nodeMoved)
                    {
                        xyz[node] = previous;
                        continue;
                    }
                    movedNodes.Add(node);
                    previousXyz.Add(previous);
                    double distance2 = 0;
                    for (int d = 0; d < 3; d++) distance2 += (point[d] - previous[d]) * (point[d] - previous[d]);
                    if (!(Math.Sqrt(distance2) > RegionMoveTolerance * meanLength)) continue;
                    moved = true;
                    for (int k = cornerStart[i]; k < cornerStart[i + 1]; k++)
                    {
                        foreach (var elementNode in elements[cornerElements[k]])
                        {
                            if (nodeIndex.TryGetValue(elementNode, out index)) changed[index] = true;
                        }
                    }
                }
                surface.Project(xyz, normals, movedNodes, previousXyz.ToArray());
                movedNodes.Clear();
                previousXyz.Clear();
                if (!moved) break;
            }
            return GetSpaceMinScaledJacobian(xyz, normals, elements, orientation);
        }
        /// <summary>
        /// Two unit vectors of the plane with the given unit normal, at right angles, that turn with it.
        /// </summary>
        private static void GetTangents(double[] normal, double[] e1, double[] e2)
        {
            int axis = 0;
            for (int d = 1; d < 3; d++)
            {
                if (Math.Abs(normal[d]) < Math.Abs(normal[axis])) axis = d;
            }
            double[] a = new double[3];
            a[axis] = 1;
            e1[0] = a[1] * normal[2] - a[2] * normal[1];
            e1[1] = a[2] * normal[0] - a[0] * normal[2];
            e1[2] = a[0] * normal[1] - a[1] * normal[0];
            double length = Math.Sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
            for (int d = 0; d < 3; d++) e1[d] /= length;
            e2[0] = normal[1] * e1[2] - normal[2] * e1[1];
            e2[1] = normal[2] * e1[0] - normal[0] * e1[2];
            e2[2] = normal[0] * e1[1] - normal[1] * e1[0];
        }
        /// <summary>
        /// Turns the two vectors of the tangent plane of a free node so that the first points the way its measure
        /// falls fastest - the condition energy of the corners it changes, or the negative of their worst scaled
        /// Jacobian - by the differences of the measure over the length h either way along the two. The directions a
        /// node is then moved in belong to the mesh alone: turned in space or numbered another way it is smoothed
        /// the same. The vectors stay as they are where the measure is not finite or does not change.
        /// </summary>
        private static void TurnTangentsDownhill(double[][] xyz, double[][] normals, int[][] elements,
                                                 int[] cornerElements, int[] cornerIndices, int start, int end,
                                                 double orientation, double[] point, double[] normal, double h,
                                                 bool byEnergy, double[] e1, double[] e2)
        {
            double[] origin = (double[])point.Clone();
            double[] slope = new double[2];
            for (int k = 0; k < 2; k++)
            {
                double[] e = k == 0 ? e1 : e2;
                double[] values = new double[2];
                for (int side = 0; side < 2; side++)
                {
                    double x = side == 0 ? h : -h;
                    for (int d = 0; d < 3; d++) point[d] = origin[d] + x * e[d];
                    double energy;
                    double minCorner = GetSpaceRegionNodeQuality(xyz, normals, elements, cornerElements, cornerIndices,
                                                                 start, end, orientation, out energy);
                    values[side] = byEnergy ? energy : -minCorner;
                }
                slope[k] = values[0] - values[1];
            }
            Array.Copy(origin, point, 3);
            double length = Math.Sqrt(slope[0] * slope[0] + slope[1] * slope[1]);
            if (!(length > 0) || double.IsInfinity(length)) return;
            double[] first = new double[3];
            for (int d = 0; d < 3; d++) first[d] = -(slope[0] * e1[d] + slope[1] * e2[d]) / length;
            Array.Copy(first, e1, 3);
            e2[0] = normal[1] * e1[2] - normal[2] * e1[1];
            e2[1] = normal[2] * e1[0] - normal[0] * e1[2];
            e2[2] = normal[0] * e1[1] - normal[1] * e1[0];
        }
        private static void SetTangentPoint(double[] point, double[] origin, double[] e1, double[] e2, double x, double y)
        {
            for (int d = 0; d < 3; d++) point[d] = origin[d] + x * e1[d] + y * e2[d];
        }
        private static double GetSpaceMeanLength(double[][] xyz, double[] point, int[] neighbours)
        {
            double length = 0;
            foreach (var neighbour in neighbours)
            {
                double dx = xyz[neighbour][0] - point[0];
                double dy = xyz[neighbour][1] - point[1];
                double dz = xyz[neighbour][2] - point[2];
                length += Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }
            return length / neighbours.Length;
        }
        private static double GetSpaceRegionNodeQuality(double[][] xyz, double[][] normals, int[][] elements,
                                                        int[] cornerElements, int[] cornerIndices, int start, int end,
                                                        double orientation, out double energy)
        {
            double min = double.MaxValue;
            energy = 0;
            for (int k = start; k < end; k++)
            {
                double cornerEnergy;
                double scaledJacobian = GetSpaceCornerQuality(xyz, normals, elements[cornerElements[k]], cornerIndices[k],
                                                              orientation, out cornerEnergy);
                if (scaledJacobian < min) min = scaledJacobian;
                energy += cornerEnergy;
            }
            return min;
        }
        /// <summary>
        /// The scaled Jacobian of corner c of a quad or a triangle in space and its condition energy, as
        /// GetCornerQuality gives them in a plane, the area of the corner taken along the normal of the surface at its
        /// node times the orientation; at a node without a normal along the normal of the element itself, of a quad
        /// that of its diagonals.
        /// </summary>
        public static double GetSpaceCornerQuality(double[][] xyz, double[][] normals, int[] element, int c,
                                                   double orientation, out double energy)
        {
            int n = element.Length;
            double[] p = xyz[element[c]];
            double[] next = xyz[element[(c + 1) % n]];
            double[] previous = xyz[element[(c + n - 1) % n]];
            double ax = next[0] - p[0];
            double ay = next[1] - p[1];
            double az = next[2] - p[2];
            double bx = previous[0] - p[0];
            double by = previous[1] - p[1];
            double bz = previous[2] - p[2];
            double cx = ay * bz - az * by;
            double cy = az * bx - ax * bz;
            double cz = ax * by - ay * bx;
            double[] normal = normals[element[c]];
            double cross;
            if (normal != null) cross = orientation * (cx * normal[0] + cy * normal[1] + cz * normal[2]);
            else if (n == 4)
            {
                double[] p0 = xyz[element[0]];
                double[] p1 = xyz[element[1]];
                double[] p2 = xyz[element[2]];
                double[] p3 = xyz[element[3]];
                double mx = (p2[1] - p0[1]) * (p3[2] - p1[2]) - (p2[2] - p0[2]) * (p3[1] - p1[1]);
                double my = (p2[2] - p0[2]) * (p3[0] - p1[0]) - (p2[0] - p0[0]) * (p3[2] - p1[2]);
                double mz = (p2[0] - p0[0]) * (p3[1] - p1[1]) - (p2[1] - p0[1]) * (p3[0] - p1[0]);
                double length = Math.Sqrt(mx * mx + my * my + mz * mz);
                cross = length > 0 ? (cx * mx + cy * my + cz * mz) / length : 0;
            }
            else cross = Math.Sqrt(cx * cx + cy * cy + cz * cz);
            double aa = ax * ax + ay * ay + az * az;
            double bb = bx * bx + by * by + bz * bz;
            double lengths = Math.Sqrt(aa * bb);
            if (!(cross > 0) || !(lengths > 0))
            {
                energy = double.PositiveInfinity;
                return lengths > 0 ? cross / lengths * (n == 3 ? TriangleScale : 1) : -1;
            }
            if (n == 3)
            {
                energy = (aa + bb - (ax * bx + ay * by + az * bz)) / (SqrtThree * cross);
                return TriangleScale * cross / lengths;
            }
            energy = (aa + bb) / (2 * cross);
            return cross / lengths;
        }
        /// <summary>
        /// The worst scaled Jacobian in space of the corners of the quads and triangles, 1 for none.
        /// </summary>
        public static double GetSpaceMinScaledJacobian(double[][] xyz, double[][] normals, int[][] elements,
                                                       double orientation)
        {
            double min = 1;
            double energy;
            foreach (var element in elements)
            {
                for (int c = 0; c < element.Length; c++)
                    min = Math.Min(min, GetSpaceCornerQuality(xyz, normals, element, c, orientation, out energy));
            }
            return min;
        }
        private static void AddUnique(List<int> list, int value)
        {
            if (!list.Contains(value)) list.Add(value);
        }
        private static double GetMeanLength(double[][] points, double[] point, int[] neighbours)
        {
            double length = 0;
            foreach (var neighbour in neighbours)
            {
                double dx = points[neighbour][0] - point[0];
                double dy = points[neighbour][1] - point[1];
                length += Math.Sqrt(dx * dx + dy * dy);
            }
            return length / neighbours.Length;
        }
        /// <summary>
        /// The worst scaled Jacobian of the corners a free node of a region changes, and their condition energy, infinite
        /// when one of them folds.
        /// </summary>
        private static double GetRegionNodeQuality(double[][] points, int[][] elements, int[] cornerElements,
                                                   int[] cornerIndices, int start, int end, double orientation,
                                                   out double energy)
        {
            double min = double.MaxValue;
            energy = 0;
            for (int k = start; k < end; k++)
            {
                double cornerEnergy;
                double scaledJacobian = GetCornerQuality(points, elements[cornerElements[k]], cornerIndices[k], orientation,
                                                         out cornerEnergy);
                if (scaledJacobian < min) min = scaledJacobian;
                energy += cornerEnergy;
            }
            return min;
        }
        /// <summary>
        /// The scaled Jacobian of corner c of a quad or a triangle and its condition energy, infinite where the corner
        /// turns against the orientation or has no area.
        /// </summary>
        private static double GetCornerQuality(double[][] points, int[] element, int c, double orientation,
                                               out double energy)
        {
            int n = element.Length;
            double[] p = points[element[c]];
            double[] next = points[element[(c + 1) % n]];
            double[] previous = points[element[(c + n - 1) % n]];
            double ax = next[0] - p[0];
            double ay = next[1] - p[1];
            double bx = previous[0] - p[0];
            double by = previous[1] - p[1];
            double cross = orientation * (ax * by - ay * bx);
            double aa = ax * ax + ay * ay;
            double bb = bx * bx + by * by;
            double lengths = Math.Sqrt(aa * bb);
            if (!(cross > 0) || !(lengths > 0))
            {
                energy = double.PositiveInfinity;
                return lengths > 0 ? cross / lengths * (n == 3 ? TriangleScale : 1) : -1;
            }
            if (n == 3)
            {
                energy = (aa + bb - (ax * bx + ay * by)) / (SqrtThree * cross);
                return TriangleScale * cross / lengths;
            }
            energy = (aa + bb) / (2 * cross);
            return cross / lengths;
        }
        /// <summary>
        /// The scaled Jacobian of corner c of a quad or a triangle, a triangle scaled so that its corner of 60 degrees
        /// counts as 1.
        /// </summary>
        public static double GetElementScaledJacobian(double[][] points, int[] element, int c, double orientation)
        {
            double energy;
            return GetCornerQuality(points, element, c, orientation, out energy);
        }
        /// <summary>
        /// The worst scaled Jacobian of the corners of the quads and triangles, 1 for none.
        /// </summary>
        public static double GetMinScaledJacobian(double[][] points, int[][] elements, double orientation)
        {
            double min = 1;
            foreach (var element in elements)
            {
                for (int c = 0; c < element.Length; c++)
                    min = Math.Min(min, GetElementScaledJacobian(points, element, c, orientation));
            }
            return min;
        }
        /// <summary>
        /// The number of quads and triangles of which a corner folds or has no area.
        /// </summary>
        public static int CountFoldedElements(double[][] points, int[][] elements, double orientation)
        {
            int numOfFolded = 0;
            foreach (var element in elements)
            {
                for (int c = 0; c < element.Length; c++)
                {
                    if (GetElementScaledJacobian(points, element, c, orientation) <= 0)
                    {
                        numOfFolded++;
                        break;
                    }
                }
            }
            return numOfFolded;
        }
        internal static double[][] GetSearchDirections()
        {
            double diagonal = Math.Sqrt(0.5);
            return new double[][]
            {
                new double[] { 1, 0 }, new double[] { -1, 0 }, new double[] { 0, 1 }, new double[] { 0, -1 },
                new double[] { diagonal, diagonal }, new double[] { -diagonal, diagonal },
                new double[] { diagonal, -diagonal }, new double[] { -diagonal, -diagonal }
            };
        }
        /// <summary>
        /// The sine of the angle at corner c of the quad, negative where the quad turns against the orientation.
        /// </summary>
        internal static double GetScaledJacobian(double[][] points, int[] quad, int c, double orientation)
        {
            double[] p = points[quad[c]];
            double[] next = points[quad[(c + 1) % 4]];
            double[] previous = points[quad[(c + 3) % 4]];
            double ax = next[0] - p[0];
            double ay = next[1] - p[1];
            double bx = previous[0] - p[0];
            double by = previous[1] - p[1];
            double lengths = Math.Sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
            if (!(lengths > 0)) return -1;
            return orientation * (ax * by - ay * bx) / lengths;
        }
        /// <summary>
        /// The negative condition energy of the corners a node changes, minus infinity when one of them folds.
        /// </summary>
        internal static double GetNodeQuality(double[][] points, int[][] quads, int[] cornerStart, int[] cornerQuads,
                                              int[] cornerIndices, int node, double orientation)
        {
            double energy = 0;
            for (int i = cornerStart[node]; i < cornerStart[node + 1]; i++)
            {
                double cornerEnergy = GetConditionEnergy(points, quads[cornerQuads[i]], cornerIndices[i], orientation);
                if (double.IsPositiveInfinity(cornerEnergy)) return double.NegativeInfinity;
                energy += cornerEnergy;
            }
            return -energy;
        }
        /// <summary>
        /// (|a|^2 + |b|^2) / (2 a x b) of corner c of the quad, infinite where the corner turns against the orientation
        /// or has no area.
        /// </summary>
        private static double GetConditionEnergy(double[][] points, int[] quad, int c, double orientation)
        {
            double[] p = points[quad[c]];
            double[] next = points[quad[(c + 1) % 4]];
            double[] previous = points[quad[(c + 3) % 4]];
            double ax = next[0] - p[0];
            double ay = next[1] - p[1];
            double bx = previous[0] - p[0];
            double by = previous[1] - p[1];
            double cross = orientation * (ax * by - ay * bx);
            if (!(cross > 0)) return double.PositiveInfinity;
            return (ax * ax + ay * ay + bx * bx + by * by) / (2 * cross);
        }
        private static double GetMinScaledJacobian(double[][] points, int[][] quads, int[] cornerStart, int[] cornerQuads,
                                                   int[] cornerIndices, int node, double orientation)
        {
            double min = double.MaxValue;
            for (int i = cornerStart[node]; i < cornerStart[node + 1]; i++)
                min = Math.Min(min, GetScaledJacobian(points, quads[cornerQuads[i]], cornerIndices[i], orientation));
            return min;
        }
        private static bool IsAnyQuadFolded(double[][] points, int[][] quads, int[] nodeQuads, int start, int end,
                                            double orientation)
        {
            for (int i = start; i < end; i++)
            {
                int[] quad = quads[nodeQuads[i]];
                for (int c = 0; c < 4; c++)
                {
                    if (GetScaledJacobian(points, quad, c, orientation) <= 0) return true;
                }
            }
            return false;
        }
        internal static double[][] CopyPoints(double[][] points)
        {
            double[][] copy = new double[points.Length][];
            for (int n = 0; n < points.Length; n++) copy[n] = new double[] { points[n][0], points[n][1] };
            return copy;
        }
    }
}
