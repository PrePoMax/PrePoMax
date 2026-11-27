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

namespace CaeMesh
{
    /// <summary>
    /// Places the holes of a submapped face in its lattice once the numbers of elements of all its edges are known, the
    /// same way in the detection and in the meshers. The offsets of the holes are not unknowns of the interval
    /// problem: the counts only make room for them (SubmapFace.GetRoomRows), and a hole is then put where it lies.
    /// Every node of a hole has a lattice position in the grid of the outer loop alone (SubmapMapping.Map, with
    /// straight lattice lines where they give a good grid), in the 2D coordinates of the nodes, the parameters of a
    /// mesher or the plane of the detection: the cell of that grid the node lies in and its place in the cell, by the
    /// inverse of the bilinear map of the cell. The target of the hole is the mean, over its nodes, of the offset
    /// that would put the node where it lies. So the hole lies between the lattice lines that pass it, also where
    /// they are not evenly spaced, whatever its shape, and exactly where they pass through its nodes. The hole is
    /// offset as is nearest to the target, and then, offset by offset in the order of
    /// their distances from the target, ties by i and then by j, until the loops of the face are valid
    /// (SubmapMapping.AreCountsValid) and the hole keeps at least the bound of cells from the outer outline and from
    /// the holes placed before it. The holes are placed one after another in loop order.
    /// A round hole is a square of a quarter of its nodes on each side, placed as any hole; where its corners go is
    /// found before, by the lattice positions of its nodes about its centre (GetRoundLoopOrder), the same way in the
    /// detection and in the meshers, which then pass its nodes in that order. Offset distances within
    /// OffsetTieTolerance are tied and ordered by i then j.
    /// </summary>
    public static class SubmapHolePlacement
    {
        // Variables                                                                                                                
        /// <summary>
        /// How far from the offset nearest to its target a hole is looked for, in lattice cells.
        /// </summary>
        public const int MaxSearchRadius = 12;
        private const double OffsetTieTolerance = 1E-9;
        private const double AngleTieTolerance = 1E-9;
        private const int MaxCellIterations = 8;
        private static readonly int[] DirectionI = new int[] { 1, 0, -1, 0 };
        private static readonly int[] DirectionJ = new int[] { 0, 1, 0, -1 };


        // Methods                                                                                                                  
        /// <summary>
        /// The lattice offsets of the first nodes of the holes, loops 1 and on, from the first node of the outer loop,
        /// loop 0. loopNodePoints holds the 2D coordinates of the nodes of every loop in loop order, node i of edge k
        /// after the nodes of the edges before it, one per element. False, with the reason, when a loop is not valid on
        /// its own or no place is found for a hole; the reason names the edges by their places in all loops one after
        /// another.
        /// </summary>
        public static bool TryPlace(int[][] loopDirections, int[][] loopCounts, double[][][] loopNodePoints, int bound,
                                    out int[][] holeOffsets, out string reason)
        {
            holeOffsets = null;
            if (!SubmapMapping.AreCountsValid(loopDirections, loopCounts, null, out reason)) return false;
            int numOfLoops = loopDirections.Length;
            holeOffsets = new int[numOfLoops - 1][];
            if (numOfLoops == 1) return true;
            //
            int[][] loopI = new int[numOfLoops][];
            int[][] loopJ = new int[numOfLoops][];
            int[] firstEdges = new int[numOfLoops];
            for (int l = 0; l < numOfLoops; l++)
            {
                if (l > 0) firstEdges[l] = firstEdges[l - 1] + loopDirections[l - 1].Length;
                GetLatticeCoordinates(loopDirections[l], loopCounts[l], out loopI[l], out loopJ[l]);
                if (loopNodePoints == null || loopNodePoints.Length != numOfLoops || loopNodePoints[l] == null ||
                    loopNodePoints[l].Length != loopI[l].Length)
                {
                    reason = "The loops of the face have other numbers of nodes than of elements.";
                    return false;
                }
            }
            double[][] outerPoints;
            SubmapGrid outerGrid;
            try
            {
                outerGrid = MapOuterLoop(loopDirections[0], loopCounts[0], loopNodePoints[0], out outerPoints);
            }
            catch (CaeException ex)
            {
                reason = ex.Message;
                holeOffsets = null;
                return false;
            }
            //
            for (int h = 1; h < numOfLoops; h++)
            {
                double targetI = 0;
                double targetJ = 0;
                bool located = true;
                for (int p = 0; p < loopI[h].Length && located; p++)
                {
                    double[] position = GetLatticePosition(outerGrid, outerPoints, loopNodePoints[h][p]);
                    located = position != null;
                    if (!located) break;
                    targetI += (position[0] - loopI[h][p]) / loopI[h].Length;
                    targetJ += (position[1] - loopJ[h][p]) / loopI[h].Length;
                }
                if (!located)
                {
                    reason = "The hole starting at Edge " + firstEdges[h] + " lies in no cell of the grid of the " +
                             "outer loop.";
                    holeOffsets = null;
                    return false;
                }
                //
                int[][] placedDirections = new int[h + 1][];
                int[][] placedCounts = new int[h + 1][];
                int[][] placedOffsets = new int[h][];
                for (int l = 0; l <= h; l++)
                {
                    placedDirections[l] = loopDirections[l];
                    placedCounts[l] = loopCounts[l];
                    if (l > 0 && l < h) placedOffsets[l - 1] = holeOffsets[l - 1];
                }
                int[] offset = null;
                foreach (var candidate in GetCandidateOffsets(targetI, targetJ))
                {
                    placedOffsets[h - 1] = candidate;
                    string candidateReason;
                    if (!SubmapMapping.AreCountsValid(placedDirections, placedCounts, placedOffsets, out candidateReason))
                        continue;
                    if (bound > 1 && !KeepsBound(loopI, loopJ, holeOffsets, h, candidate, bound)) continue;
                    offset = candidate;
                    break;
                }
                if (offset == null)
                {
                    reason = "The lattice of the face has no place for the hole starting at Edge " + firstEdges[h] +
                             " that keeps " + (bound == 1 ? "a cell" : bound + " cells") + " between it and the other " +
                             "loops.";
                    holeOffsets = null;
                    return false;
                }
                holeOffsets[h - 1] = offset;
            }
            return true;
        }
        /// <summary>
        /// The grid of the outer loop alone on its nodes, one per element in loop order, and the 2D coordinates of all
        /// its nodes, the nodes of the loop first.
        /// </summary>
        private static SubmapGrid MapOuterLoop(int[] directions, int[] counts, double[][] nodePoints,
                                               out double[][] points)
        {
            int numOfNodes = nodePoints.Length;
            double[][][] edgeUv = new double[directions.Length][][];
            double[][][] edgeXyz = new double[directions.Length][][];
            int first = 0;
            for (int k = 0; k < directions.Length; k++)
            {
                edgeUv[k] = new double[counts[k] + 1][];
                edgeXyz[k] = new double[counts[k] + 1][];
                for (int s = 0; s <= counts[k]; s++)
                {
                    double[] point = nodePoints[(first + s) % numOfNodes];
                    edgeUv[k][s] = point;
                    edgeXyz[k][s] = new double[] { point[0], point[1], 0 };
                }
                first += counts[k];
            }
            SubmapGrid grid = SubmapMapping.Map(new int[][] { directions }, new double[][][][] { edgeUv },
                                                new double[][][][] { edgeXyz }, new int[0][], null, false, true, null);
            points = nodePoints.Concat(grid.InnerUv).ToArray();
            return grid;
        }
        /// <summary>
        /// The lattice coordinates of the point from the first node of the grid: those of the cell whose bilinear map
        /// puts the point nearest to its middle, among the cells no farther from the point than they are large, and
        /// the place of the point in it. Null when no cell gives a place.
        /// </summary>
        private static double[] GetLatticePosition(SubmapGrid grid, double[][] points, double[] point)
        {
            int[][] logical = grid.LogicalCoordinates;
            double[] position = null;
            double least = double.PositiveInfinity;
            foreach (var quad in grid.Quads)
            {
                double[] a = points[quad[0]];
                double[] b = points[quad[1]];
                double[] c = points[quad[2]];
                double[] d = points[quad[3]];
                double minX = Math.Min(Math.Min(a[0], b[0]), Math.Min(c[0], d[0]));
                double maxX = Math.Max(Math.Max(a[0], b[0]), Math.Max(c[0], d[0]));
                double minY = Math.Min(Math.Min(a[1], b[1]), Math.Min(c[1], d[1]));
                double maxY = Math.Max(Math.Max(a[1], b[1]), Math.Max(c[1], d[1]));
                double margin = Math.Max(maxX - minX, maxY - minY);
                if (point[0] < minX - margin || point[0] > maxX + margin || point[1] < minY - margin ||
                    point[1] > maxY + margin) continue;
                double s = 0.5;
                double t = 0.5;
                for (int iteration = 0; iteration < MaxCellIterations; iteration++)
                {
                    double fx = (1 - s) * (1 - t) * a[0] + s * (1 - t) * b[0] + s * t * c[0] + (1 - s) * t * d[0] -
                                point[0];
                    double fy = (1 - s) * (1 - t) * a[1] + s * (1 - t) * b[1] + s * t * c[1] + (1 - s) * t * d[1] -
                                point[1];
                    double xs = (1 - t) * (b[0] - a[0]) + t * (c[0] - d[0]);
                    double ys = (1 - t) * (b[1] - a[1]) + t * (c[1] - d[1]);
                    double xt = (1 - s) * (d[0] - a[0]) + s * (c[0] - b[0]);
                    double yt = (1 - s) * (d[1] - a[1]) + s * (c[1] - b[1]);
                    double determinant = xs * yt - xt * ys;
                    s -= (fx * yt - fy * xt) / determinant;
                    t -= (fy * xs - fx * ys) / determinant;
                }
                double distance = Math.Max(Math.Abs(s - 0.5), Math.Abs(t - 0.5));
                if (!(distance < least)) continue;
                least = distance;
                position = new double[] { logical[quad[0]][0] - logical[0][0] + s,
                                          logical[quad[0]][1] - logical[0][1] + t };
            }
            return position;
        }
        /// <summary>
        /// Whether the holes of a face find a place with the given counts of its edges, as the balance of the counts
        /// sees the face before it is meshed: its nodes laid in equal steps along straight logical edges between its
        /// vertices, which loopVertexPoints holds as SubmapFace.GetLogicalLoopPoints takes them, and the nodes of a round
        /// hole evenly spread over its circle and started by GetRoundLoopOrder. True for a face without holes; the reason
        /// names the edges by their places in all loops one after another.
        /// </summary>
        public static bool TryPlace(SubmapFace face, Dictionary<int, int> edgeIdCount, double[][][] loopVertexPoints,
                                    int bound, out string reason)
        {
            reason = null;
            if (!face.HasHoles) return true;
            int[][] counts = face.GetLoopLogicalEdgeCounts(edgeIdCount);
            double[][][] logicalPoints = face.GetLogicalLoopPoints(loopVertexPoints);
            if (logicalPoints == null)
            {
                reason = "The vertices of the face give no frame for its holes.";
                return false;
            }
            double[][][] nodePoints = new double[counts.Length][][];
            SubmapFace[] loops = face.Loops;
            int[][] directions = face.LoopLogicalEdgeDirections;
            for (int l = 0; l < counts.Length; l++)
            {
                if (!loops[l].IsRound)
                {
                    nodePoints[l] = GetStraightNodePoints(logicalPoints[l], counts[l]);
                    continue;
                }
                double[][] circlePoints = GetCircleNodePoints(logicalPoints[l], 4 * counts[l][0]);
                int[] order = GetRoundLoopOrder(directions[0], counts[0], nodePoints[0], circlePoints);
                nodePoints[l] = new double[order.Length][];
                for (int p = 0; p < order.Length; p++) nodePoints[l][p] = circlePoints[order[p]];
            }
            int[][] holeOffsets;
            return TryPlace(face.LoopLogicalEdgeDirections, counts, nodePoints, bound, out holeOffsets, out reason);
        }
        /// <summary>
        /// The 2D coordinates of the nodes of a loop whose logical edge k runs straight from logicalVertexPoints[k] to
        /// the next one in counts[k] equal steps, one point per node in loop order.
        /// </summary>
        public static double[][] GetStraightNodePoints(double[][] logicalVertexPoints, int[] counts)
        {
            int n = counts.Length;
            List<double[]> points = new List<double[]>();
            for (int k = 0; k < n; k++)
            {
                double[] a = logicalVertexPoints[k];
                double[] b = logicalVertexPoints[(k + 1) % n];
                for (int s = 0; s < counts[k]; s++)
                {
                    double f = (double)s / counts[k];
                    points.Add(new double[] { a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f });
                }
            }
            return points.ToArray();
        }
        /// <summary>
        /// The 2D coordinates of the given number of nodes evenly spread over the circle of a round hole, clockwise from
        /// where it reaches furthest along +i; roundLogicalPoints are the logical points of the hole that
        /// SubmapFace.GetLogicalLoopPoints gives, whose last four span the circle.
        /// </summary>
        private static double[][] GetCircleNodePoints(double[][] roundLogicalPoints, int numOfNodes)
        {
            double[] east = roundLogicalPoints[4];
            double[] west = roundLogicalPoints[6];
            double[] centre = new double[] { 0.5 * (east[0] + west[0]), 0.5 * (east[1] + west[1]) };
            double radius = 0.5 * (east[0] - west[0]);
            double[][] points = new double[numOfNodes][];
            for (int p = 0; p < numOfNodes; p++)
            {
                double angle = -2 * Math.PI * p / numOfNodes;
                points[p] = new double[] { centre[0] + radius * Math.Cos(angle), centre[1] + radius * Math.Sin(angle) };
            }
            return points;
        }
        /// <summary>
        /// The order in which the nodes of a round hole make its logical loop, node p of the logical loop being node
        /// order[p] of loopNodePoints, which holds the 2D coordinates of its nodes as they follow each other around it,
        /// in either sense. The outer loop gives the lattice: its logical edges run in outerDirections with outerCounts
        /// elements, through the nodes outerNodePoints, one per element. Every node of the hole and its centre, the mean
        /// of its nodes, are given their lattice positions by the mean value coordinates interpolation of the lattice
        /// coordinates of the outer nodes, as the target of a hole is; the logical loop starts at the +i +j corner of
        /// the square of the hole and runs the way that is clockwise in the lattice, so that its four logical edges of
        /// a quarter of the nodes each run along -j, -i, +j and +i. The corners are the four nodes, a quarter of the
        /// nodes apart, whose angles about the centre are nearest to the diagonals of the lattice, by the sum of their
        /// squares: no corner is told from the others, so the mirror image of a hole gets the mirror image of its
        /// square. A hole that is its own mirror image with a node in the plane and an odd number of elements on a side
        /// of its square has two such squares, turned either way by the same angle: the one is taken that is turned
        /// the way the place of the hole in the lattice box of the outer loop tells, which is the other way in a
        /// mirror image of the place, and by the smaller signed angle from 45° where the place tells nothing, on a
        /// middle line or a diagonal of the box. The sense is told by the angles of the nodes, not by their order: the
        /// sum of the turns from node to node about the centre.
        /// </summary>
        public static int[] GetRoundLoopOrder(int[] outerDirections, int[] outerCounts, double[][] outerNodePoints,
                                              double[][] loopNodePoints)
        {
            int[] outerI;
            int[] outerJ;
            GetLatticeCoordinates(outerDirections, outerCounts, out outerI, out outerJ);
            double[][] outerValues = new double[outerI.Length][];
            for (int p = 0; p < outerValues.Length; p++) outerValues[p] = new double[] { outerI[p], outerJ[p] };
            int n = loopNodePoints.Length;
            double[] centre = new double[2];
            foreach (var point in loopNodePoints)
            {
                centre[0] += point[0];
                centre[1] += point[1];
            }
            centre[0] /= n;
            centre[1] /= n;
            double[] target = InterpolateMeanValue(outerNodePoints, outerValues, centre);
            double[] angles = new double[n];
            for (int p = 0; p < n; p++)
            {
                double[] position = InterpolateMeanValue(outerNodePoints, outerValues, loopNodePoints[p]);
                angles[p] = Math.Atan2(position[1] - target[1], position[0] - target[0]);
            }
            double sumOfTurns = 0;
            for (int p = 0; p < n; p++) sumOfTurns += GetWrappedAngle(angles[(p + 1) % n] - angles[p]);
            int step = sumOfTurns > 0 ? -1 : 1;
            // How far the corners of the square that starts at each node lie from the diagonals of the lattice
            int numOfCorners = n % 4 == 0 ? 4 : 1;
            double[] deviations = new double[n];
            double[] sums = new double[n];
            double least = double.PositiveInfinity;
            for (int p = 0; p < n; p++)
            {
                deviations[p] = GetWrappedAngle(angles[p] - 0.25 * Math.PI);
                for (int k = 0; k < numOfCorners; k++)
                {
                    int corner = ((p + step * k * (n / 4)) % n + n) % n;
                    double deviation = GetWrappedAngle(angles[corner] - (0.25 - 0.5 * k) * Math.PI);
                    sums[p] += deviation * deviation;
                }
                least = Math.Min(least, sums[p]);
            }
            // The way a tied square is turned, by the place of the hole from the middle of the lattice box
            double placeI = target[0] - 0.5 * (outerI.Min() + outerI.Max());
            double placeJ = target[1] - 0.5 * (outerJ.Min() + outerJ.Max());
            double turn = placeI * placeJ * (placeI * placeI - placeJ * placeJ);
            bool larger = turn < -AngleTieTolerance * Math.Max(1, Math.Pow(Math.Abs(placeI) + Math.Abs(placeJ), 4));
            int start = -1;
            for (int p = 0; p < n; p++)
            {
                if (sums[p] > least + AngleTieTolerance) continue;
                if (start < 0 || (larger ? deviations[p] > deviations[start] : deviations[p] < deviations[start]))
                    start = p;
            }
            int[] order = new int[n];
            for (int p = 0; p < n; p++) order[p] = ((start + step * p) % n + n) % n;
            return order;
        }
        /// <summary>
        /// The angle brought into -π to π.
        /// </summary>
        private static double GetWrappedAngle(double angle)
        {
            while (angle > Math.PI) angle -= 2 * Math.PI;
            while (angle <= -Math.PI) angle += 2 * Math.PI;
            return angle;
        }
        /// <summary>
        /// The mean value coordinates interpolation, at the point, of the values given at the vertices of a closed
        /// polygon. At a vertex its value, and on an edge the linear interpolation between its two vertices.
        /// </summary>
        public static double[] InterpolateMeanValue(double[][] polygon, double[][] values, double[] point)
        {
            int n = polygon.Length;
            int dimension = values[0].Length;
            double[][] d = new double[n][];
            double[] r = new double[n];
            double size = 0;
            for (int i = 0; i < n; i++)
            {
                d[i] = new double[] { polygon[i][0] - point[0], polygon[i][1] - point[1] };
                r[i] = Math.Sqrt(d[i][0] * d[i][0] + d[i][1] * d[i][1]);
                size = Math.Max(size, r[i]);
            }
            double tolerance = 1E-12 * size;
            for (int i = 0; i < n; i++)
            {
                if (r[i] <= tolerance) return (double[])values[i].Clone();
            }
            // The tangent of half the angle at the point from vertex i to vertex i + 1
            double[] tangents = new double[n];
            for (int i = 0; i < n; i++)
            {
                int next = (i + 1) % n;
                double cross = d[i][0] * d[next][1] - d[i][1] * d[next][0];
                double dot = d[i][0] * d[next][0] + d[i][1] * d[next][1];
                if (Math.Abs(cross) <= tolerance * Math.Max(r[i], r[next]) && dot < 0)
                {
                    double f = r[i] / (r[i] + r[next]);
                    double[] onEdge = new double[dimension];
                    for (int c = 0; c < dimension; c++) onEdge[c] = values[i][c] + (values[next][c] - values[i][c]) * f;
                    return onEdge;
                }
                double numerator = r[i] * r[next] - dot;
                tangents[i] = numerator <= tolerance * Math.Max(r[i], r[next]) ? 0 : numerator / cross;
            }
            double[] sum = new double[dimension];
            double sumOfWeights = 0;
            for (int i = 0; i < n; i++)
            {
                double weight = (tangents[(i + n - 1) % n] + tangents[i]) / r[i];
                for (int c = 0; c < dimension; c++) sum[c] += weight * values[i][c];
                sumOfWeights += weight;
            }
            for (int c = 0; c < dimension; c++) sum[c] /= sumOfWeights;
            return sum;
        }
        /// <summary>
        /// The lattice coordinates of the nodes of a loop from its first node, one per element in loop order.
        /// </summary>
        private static void GetLatticeCoordinates(int[] directions, int[] counts, out int[] nodeI, out int[] nodeJ)
        {
            int numOfNodes = 0;
            foreach (var count in counts) numOfNodes += count;
            nodeI = new int[numOfNodes];
            nodeJ = new int[numOfNodes];
            int i = 0;
            int j = 0;
            int p = 0;
            for (int k = 0; k < directions.Length; k++)
            {
                for (int s = 0; s < counts[k]; s++)
                {
                    nodeI[p] = i;
                    nodeJ[p] = j;
                    p++;
                    i += DirectionI[directions[k]];
                    j += DirectionJ[directions[k]];
                }
            }
        }
        /// <summary>
        /// The integer offsets within MaxSearchRadius of the target, by their distances from it, ties by i and then j.
        /// </summary>
        private static List<int[]> GetCandidateOffsets(double targetI, double targetJ)
        {
            int centreI = GetNearestOffset(targetI);
            int centreJ = GetNearestOffset(targetJ);
            List<int[]> candidates = new List<int[]>();
            List<double> distances = new List<double>();
            for (int di = -MaxSearchRadius; di <= MaxSearchRadius; di++)
            {
                for (int dj = -MaxSearchRadius; dj <= MaxSearchRadius; dj++)
                {
                    if (di * di + dj * dj > MaxSearchRadius * MaxSearchRadius) continue;
                    int i = centreI + di;
                    int j = centreJ + dj;
                    candidates.Add(new int[] { i, j });
                    distances.Add(Math.Sqrt((i - targetI) * (i - targetI) + (j - targetJ) * (j - targetJ)));
                }
            }
            int[] order = new int[candidates.Count];
            for (int c = 0; c < order.Length; c++) order[c] = c;
            Array.Sort(order, (a, b) => distances[a].CompareTo(distances[b]));
            var byOffset = Comparer<int>.Create((a, b) =>
            {
                int compareI = candidates[a][0].CompareTo(candidates[b][0]);
                return compareI != 0 ? compareI : candidates[a][1].CompareTo(candidates[b][1]);
            });
            for (int start = 0; start < order.Length;)
            {
                int end = start + 1;
                double lastDistance = distances[order[start]] + OffsetTieTolerance;
                while (end < order.Length && distances[order[end]] <= lastDistance) end++;
                Array.Sort(order, start, end - start, byOffset);
                start = end;
            }
            List<int[]> sorted = new List<int[]>();
            foreach (var c in order) sorted.Add(candidates[c]);
            return sorted;
        }
        private static int GetNearestOffset(double target)
        {
            int lower = (int)Math.Floor(target);
            return target - lower > 0.5 + OffsetTieTolerance ? lower + 1 : lower;
        }
        /// <summary>
        /// Whether hole h at the offset keeps at least the bound of cells, by the lattice distance along i or j,
        /// whichever is larger, between its nodes and those of the outer loop and of the holes placed before it.
        /// </summary>
        private static bool KeepsBound(int[][] loopI, int[][] loopJ, int[][] holeOffsets, int h, int[] offset, int bound)
        {
            for (int l = 0; l < h; l++)
            {
                int offsetI = l == 0 ? 0 : holeOffsets[l - 1][0];
                int offsetJ = l == 0 ? 0 : holeOffsets[l - 1][1];
                for (int p = 0; p < loopI[h].Length; p++)
                {
                    int i = loopI[h][p] + offset[0];
                    int j = loopJ[h][p] + offset[1];
                    for (int q = 0; q < loopI[l].Length; q++)
                    {
                        int distance = Math.Max(Math.Abs(loopI[l][q] + offsetI - i), Math.Abs(loopJ[l][q] + offsetJ - j));
                        if (distance < bound) return false;
                    }
                }
            }
            return true;
        }
    }
}
