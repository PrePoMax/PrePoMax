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
    /// Pairs boundary nodes along a submap lattice. Cornered holes supply nodes and interrupt lines of sight.
    /// Lines pass through round holes. With the round holes of a plane face placed, the lattice lines that end on a
    /// round hole get targets: the even places, without the first and the last, over the hole, the nodes its grid is
    /// laid on (those of the circle on the outer side of its ring, or its own without a ring) and a clearance of
    /// RoundHoleClearance of the distance between two of them. The lattice is then centred on the hole, and the
    /// lines beside the hole, which the evening places, pass it clear.
    /// </summary>
    internal static class GmshShellSubmapSpacing
    {
        // Variables                                                                                                   
        private static readonly int[] _directionI = { 1, 0, -1, 0 };
        private static readonly int[] _directionJ = { 0, 1, 0, -1 };
        private const double RoundHoleClearance = 0.5;


        // Methods                                                                                                     
        public static List<(GmshShellSpacingNode First, GmshShellSpacingNode Second)>
            GetPairs(SubmapFace face, Dictionary<int, GmshShellSpacingEdge> edges, bool withRoundHoles,
                     out List<(GmshShellSpacingNode Node, double Target)> targets)
        {
            targets = new List<(GmshShellSpacingNode Node, double Target)>();
            if (withRoundHoles && GmshCommon.Gmsh.Model.GetType(2, face.FaceId) != "Plane") return null;
            SubmapFace[] loops = face.Loops;
            Dictionary<int, int> edgeCounts = new Dictionary<int, int>();
            int[][] counts = new int[loops.Length][];
            int[][] directions = face.LoopLogicalEdgeDirections;
            for (int l = 0; l < loops.Length; l++)
            {
                SubmapFace loop = loops[l];
                if (loop.IsRound) continue;
                for (int k = 0; k < loop.EdgeIds.Length; k++)
                {
                    int id = loop.EdgeIds[k];
                    GmshShellSpacingEdge edge = edges[id];
                    if (edge == null || (edge.NumOfNodes - 1) % loop.GetEdgeMultiple(k) != 0) return null;
                    edgeCounts.Add(id, edge.NumOfNodes - 1);
                }
                counts[l] = loop.GetLogicalEdgeCounts(edgeCounts);
            }
            double[] lengths = face.EdgeIds.Select(id => edges[id].Length).ToArray();
            double[][] outer = SubmapFaceDetector.GetIdealOutline(directions[0], face.GetLogicalEdgeLengths(lengths));
            if (outer == null) return null;
            //
            var loopNodes = new (GmshShellSpacingEdge Edge, int Index, int Step)[loops.Length][];
            double[][][] outlines = new double[loops.Length][][];
            outlines[0] = outer;
            for (int l = 0; l < loops.Length; l++)
            {
                if (loops[l].IsRound) continue;
                loopNodes[l] = ReadLoop(loops[l], edges, counts[l].Sum());
                if (loopNodes[l] == null || counts[l].Any(count => count <= 0)) return null;
            }
            int[][] offsets = null;
            double[][][] ringUv = new double[loops.Length][][];
            double[][][] ringPoints = new double[loops.Length][][];
            if (withRoundHoles ? loops.Length > 1 : counts.Skip(1).Any(count => count != null))
            {
                double[][][] uv;
                if (loops.Any(loop => loop.IsRound))
                {
                    if (!GmshSubmapMesher.TryPlaceHoles(face, out _, out double[][][][] loopXyz,
                                                        out double[][][][] loopUv, out uv, out offsets, out _))
                        return null;
                    double[][][][] gridUv = loopUv;
                    if (withRoundHoles)
                        GmshSubmapMesher.GetRingLoops(face.FaceId, loops, loopXyz, loopUv, out _, out gridUv);
                    for (int l = 1; l < loops.Length; l++)
                    {
                        if (!loops[l].IsRound || !withRoundHoles) continue;
                        counts[l] = gridUv[l].Select(side => side.Length - 1).ToArray();
                        ringUv[l] = gridUv[l].SelectMany(side => side.Take(side.Length - 1)).ToArray();
                    }
                }
                else
                {
                    uv = new double[loops.Length][][];
                    try
                    {
                        for (int l = 0; l < loops.Length; l++)
                            uv[l] = GmshThreeBlockMesher.GetParametricCoordinates(face.FaceId,
                                loopNodes[l].Select(node => node.Edge.GetPoint(node.Index)).ToList());
                    }
                    catch (CaeException)
                    {
                        return null;
                    }
                    if (!SubmapHolePlacement.TryPlace(directions, counts, uv, 1, out offsets, out _)) return null;
                }
                if (!MapHoleOutlines(uv, counts, outlines, ringUv, ringPoints)) return null;
            }
            //
            var nodes = new List<(int I, int J, GmshShellSpacingNode X, GmshShellSpacingNode Y)>();
            var roundNodes = new List<(int I, int J, double X, double Y, int Loop)>();
            var positions = new HashSet<(int I, int J)>();
            int[][] cornerI = new int[loops.Length][];
            int[][] cornerJ = new int[loops.Length][];
            for (int l = 0; l < loops.Length; l++)
            {
                if (counts[l] == null) continue;
                int i = l == 0 ? 0 : offsets[l - 1][0];
                int j = l == 0 ? 0 : offsets[l - 1][1];
                int numOfSegments = counts[l].Length;
                cornerI[l] = new int[numOfSegments + 1];
                cornerJ[l] = new int[numOfSegments + 1];
                int firstNode = 0;
                for (int k = 0; k < numOfSegments; k++)
                {
                    cornerI[l][k] = i;
                    cornerJ[l][k] = j;
                    int di = _directionI[directions[l][k]];
                    int dj = _directionJ[directions[l][k]];
                    int count = counts[l][k];
                    if (loops[l].IsRound)
                    {
                        for (int n = 0; n < count; n++, firstNode++, i += di, j += dj)
                        {
                            if (!positions.Add((i, j))) return null;
                            roundNodes.Add((i, j, ringPoints[l][firstNode][0], ringPoints[l][firstNode][1], l));
                        }
                        continue;
                    }
                    var start = loopNodes[l][firstNode];
                    double from = start.Edge.GetFraction(start.Index);
                    double span = start.Edge.GetFraction(start.Index + start.Step * count) - from;
                    if (span == 0) return null;
                    double[] a = outlines[l][k];
                    double[] b = outlines[l][(k + 1) % numOfSegments];
                    double sx = (b[0] - a[0]) / span;
                    double sy = (b[1] - a[1]) / span;
                    for (int n = 0; n < count; n++, firstNode++, i += di, j += dj)
                    {
                        if (!positions.Add((i, j))) return null;
                        var node = loopNodes[l][firstNode];
                        nodes.Add((i, j, new GmshShellSpacingNode(node.Edge, node.Index, a[0] - from * sx, sx),
                                         new GmshShellSpacingNode(node.Edge, node.Index, a[1] - from * sy, sy)));
                    }
                }
                cornerI[l][numOfSegments] = i;
                cornerJ[l][numOfSegments] = j;
                if (i != cornerI[l][0] || j != cornerJ[l][0]) return null;
                if (!loops[l].IsRound) continue;
                cornerI[l] = null;
                cornerJ[l] = null;
            }
            //
            var pairs = new List<(GmshShellSpacingNode First, GmshShellSpacingNode Second)>();
            foreach (var row in nodes.GroupBy(node => node.J))
            {
                var ordered = row.OrderBy(node => node.I).ToArray();
                for (int n = 1; n < ordered.Length; n++)
                {
                    double middle = 0.5 * (ordered[n - 1].I + ordered[n].I);
                    if (IsInterior(cornerI, cornerJ, middle, row.Key))
                        pairs.Add((ordered[n - 1].Y, ordered[n].Y));
                }
            }
            foreach (var column in nodes.GroupBy(node => node.I))
            {
                var ordered = column.OrderBy(node => node.J).ToArray();
                for (int n = 1; n < ordered.Length; n++)
                {
                    double middle = 0.5 * (ordered[n - 1].J + ordered[n].J);
                    if (IsInterior(cornerJ, cornerI, middle, column.Key))
                        pairs.Add((ordered[n - 1].X, ordered[n].X));
                }
            }
            // The lattice lines that end on a round hole aim at even places over the hole, the nodes its grid is laid
            // on and a clearance around them
            for (int l = 1; l < loops.Length; l++)
            {
                var round = roundNodes.Where(node => node.Loop == l).ToArray();
                if (round.Length == 0) continue;
                double x = round.Average(node => node.X);
                double y = round.Average(node => node.Y);
                double half = 0;
                for (int n = 0; n < round.Length; n++)
                {
                    var next = round[(n + 1) % round.Length];
                    half += (Math.Sqrt(Math.Pow(round[n].X - x, 2) + Math.Pow(round[n].Y - y, 2)) +
                             RoundHoleClearance * Math.Sqrt(Math.Pow(next.X - round[n].X, 2) +
                                                            Math.Pow(next.Y - round[n].Y, 2))) / round.Length;
                }
                int minI = round.Min(node => node.I);
                int maxI = round.Max(node => node.I);
                int minJ = round.Min(node => node.J);
                int maxJ = round.Max(node => node.J);
                for (int i = minI; i <= maxI; i++)
                    AddTargets(targets, nodes.Where(node => node.I == i).Select(node => (node.J, node.X)), minJ, maxJ,
                               x + half * (2.0 * (i - minI + 1) / (maxI - minI + 2) - 1));
                for (int j = minJ; j <= maxJ; j++)
                    AddTargets(targets, nodes.Where(node => node.J == j).Select(node => (node.I, node.Y)), minI, maxI,
                               y + half * (2.0 * (j - minJ + 1) / (maxJ - minJ + 2) - 1));
            }
            return pairs;
        }
        /// <summary>
        /// Aims the nodes of one lattice line that lie along a round hole, from the last node before the hole to the
        /// first after it.
        /// </summary>
        private static void AddTargets(List<(GmshShellSpacingNode Node, double Target)> targets,
                                       IEnumerable<(int Along, GmshShellSpacingNode Node)> line, int from, int to,
                                       double target)
        {
            var ordered = line.OrderBy(node => node.Along).ToArray();
            int first = Array.FindLastIndex(ordered, node => node.Along <= from);
            int last = Array.FindIndex(ordered, node => node.Along >= to);
            if (first < 0 || last < 0) return;
            for (int n = first; n <= last; n++) targets.Add((ordered[n].Node, target));
        }
        private static (GmshShellSpacingEdge Edge, int Index, int Step)[]
            ReadLoop(SubmapFace loop, Dictionary<int, GmshShellSpacingEdge> edges, int count)
        {
            var nodes = new (GmshShellSpacingEdge Edge, int Index, int Step)[count];
            int node = 0;
            for (int e = 0; e < loop.EdgeIds.Length; e++)
            {
                GmshShellSpacingEdge edge = edges[loop.EdgeIds[e]];
                bool reversed = loop.VertexNodeIds[e] == edge.EndVertexId;
                if (!reversed && loop.VertexNodeIds[e] != edge.StartVertexId) return null;
                if (loop.VertexNodeIds[(e + 1) % loop.EdgeIds.Length] !=
                    (reversed ? edge.StartVertexId : edge.EndVertexId)) return null;
                int last = edge.NumOfNodes - 1;
                for (int k = 0; k < last; k++)
                    nodes[node++] = (edge, reversed ? last - k : k, reversed ? -1 : 1);
            }
            return nodes;
        }
        private static bool MapHoleOutlines(double[][][] uv, int[][] counts, double[][][] outlines, double[][][] ringUv,
                                            double[][][] ringPoints)
        {
            double[][] outer = outlines[0];
            int count = outer.Length;
            double u0 = 0;
            double v0 = 0;
            double x0 = 0;
            double y0 = 0;
            for (int k = 0, n = 0; k < count; n += counts[0][k], k++)
            {
                u0 += uv[0][n][0];
                v0 += uv[0][n][1];
                x0 += outer[k][0];
                y0 += outer[k][1];
            }
            u0 /= count;
            v0 /= count;
            x0 /= count;
            y0 /= count;
            double uu = 0;
            double uvSum = 0;
            double vv = 0;
            double ux = 0;
            double vx = 0;
            double uy = 0;
            double vy = 0;
            for (int k = 0, n = 0; k < count; n += counts[0][k], k++)
            {
                double u = uv[0][n][0] - u0;
                double v = uv[0][n][1] - v0;
                double x = outer[k][0] - x0;
                double y = outer[k][1] - y0;
                uu += u * u;
                uvSum += u * v;
                vv += v * v;
                ux += u * x;
                vx += v * x;
                uy += u * y;
                vy += v * y;
            }
            double determinant = uu * vv - uvSum * uvSum;
            if (!(determinant > 1E-12 * uu * vv)) return false;
            double xu = (vv * ux - uvSum * vx) / determinant;
            double xv = (uu * vx - uvSum * ux) / determinant;
            double yu = (vv * uy - uvSum * vy) / determinant;
            double yv = (uu * vy - uvSum * uy) / determinant;
            for (int l = 1; l < outlines.Length; l++)
            {
                if (counts[l] == null) continue;
                outlines[l] = new double[counts[l].Length][];
                for (int k = 0, n = 0; k < counts[l].Length; n += counts[l][k], k++)
                {
                    double u = uv[l][n][0] - u0;
                    double v = uv[l][n][1] - v0;
                    outlines[l][k] = new double[] { x0 + xu * u + xv * v, y0 + yu * u + yv * v };
                }
                if (ringUv[l] == null) continue;
                ringPoints[l] = new double[ringUv[l].Length][];
                for (int n = 0; n < ringUv[l].Length; n++)
                {
                    double u = ringUv[l][n][0] - u0;
                    double v = ringUv[l][n][1] - v0;
                    ringPoints[l][n] = new double[] { x0 + xu * u + xv * v, y0 + yu * u + yv * v };
                }
            }
            return true;
        }
        // Both sides of a lattice line must be inside, including where it passes through a reentrant corner.
        private static bool IsInterior(int[][] along, int[][] across, double midpoint, int line)
        {
            bool below = false;
            bool above = false;
            for (int l = 0; l < along.Length; l++)
            {
                if (along[l] == null) continue;
                for (int k = 0; k < along[l].Length - 1; k++)
                {
                    if (along[l][k] != along[l][k + 1] || along[l][k] <= midpoint) continue;
                    int low = Math.Min(across[l][k], across[l][k + 1]);
                    int high = Math.Max(across[l][k], across[l][k + 1]);
                    if (low < line && line <= high) below = !below;
                    if (low <= line && line < high) above = !above;
                }
            }
            return below && above;
        }
    }
}
