// PrePoMax - Copyright (C) 2016-2026 Matej Borovinšek
//
// Licensed under the terms defined in the LICENSE file located in the root directory of this source code.
//
// Source code: https://github.com/PrePoMax/PrePoMax
//
// Author: Matej Borovinšek
// Contributors:

using System;
using CaeMesh.Sweep;

namespace CaeMesh
{
    /// <summary>
    /// The nodes of one meshed edge by their fractions of its length, for GmshShellSideSpacing. The sides the edge
    /// belongs to add where they would have each node, and a step moves the nodes towards the weighted mean of that,
    /// kept in order and no closer than a third of the even spacing. Another step evens the nodes out among the
    /// intervals the edge was meshed with. An edge that is not movable only answers.
    /// </summary>
    internal class GmshShellSpacingEdge
    {
        // Variables                                                                                                    
        private const double MinimumIntervalFactor = 1.0 / 3;
        //
        private readonly SweepSideNodes _nodes;
        private readonly double[] _fractions;
        private readonly double[] _errors;
        private readonly double[] _weights;
        private readonly double[] _levels;
        private readonly int[] _blockSizes;
        private readonly double[] _indices;
        private readonly double _minimumInterval;
        private bool _hasRequests;
        private bool _movable;
        private double[][] _proposedPoints;
        private double[] _proposedParameters;


        // Properties                                                                                                   
        public double Length { get { return _nodes.Length; } }
        public int StartVertexId { get { return _nodes.StartVertexId; } }
        public int EndVertexId { get { return _nodes.EndVertexId; } }
        public int NumOfNodes { get { return _fractions.Length; } }
        public bool HasProposedNodes { get { return _proposedPoints != null; } }


        // Constructors                                                                                                 
        public GmshShellSpacingEdge(SweepSideNodes nodes)
        {
            _nodes = nodes;
            _fractions = new double[nodes.NumOfNodes];
            for (int k = 0; k < _fractions.Length; k++) _fractions[k] = nodes.GetFraction(k);
            _minimumInterval = MinimumIntervalFactor / (_fractions.Length - 1);
            if (_fractions.Length <= 2) return;
            _errors = new double[_fractions.Length];
            _weights = new double[_fractions.Length];
            _levels = new double[_fractions.Length];
            _blockSizes = new int[_fractions.Length];
            _indices = new double[_fractions.Length];
        }


        // Methods                                                                                                      
        public double GetFraction(int node)
        {
            return _fractions[node];
        }
        public void SetMovable(bool movable)
        {
            _movable = movable && _errors != null;
        }
        public double[] GetPoint(int node)
        {
            return _proposedPoints?[node] ?? _nodes.GetShellPoint(node);
        }
        public double[] EvaluatePoint(double fraction)
        {
            return _nodes.EvaluateShellPoint(fraction, out _);
        }
        public void AddError(int node, double scale, double error)
        {
            if (!_movable || scale == 0 || node == 0 || node == _fractions.Length - 1) return;
            _errors[node] += scale * error;
            _weights[node] += scale * scale;
            _hasRequests = true;
        }
        public double Relax(double step = 0.5)
        {
            if (!_hasRequests) return 0;
            _hasRequests = false;
            //
            int last = _fractions.Length - 1;
            int blocks = 0;
            for (int k = 1; k < last; k++)
            {
                double correction = _weights[k] > 0 ? step * _errors[k] / _weights[k] : 0;
                blocks = AddLevel(blocks, k, _fractions[k] + correction);
                _errors[k] = 0;
                _weights[k] = 0;
            }
            return SetLevels(blocks);
        }
        /// <summary>
        /// Moves every node of a movable edge half the way to the middle of its two neighbours as they lie, counted in
        /// the intervals the edge was meshed with, so that the nodes keep the grading of the edge where nothing else
        /// places them. The step is the same from either end of the edge.
        /// </summary>
        public void Even()
        {
            if (!_movable) return;
            int last = _fractions.Length - 1;
            for (int k = 0; k <= last; k++) _indices[k] = GetIndex(_fractions[k]);
            int blocks = 0;
            int interval = 0;
            for (int k = 1; k < last; k++)
            {
                double index = 0.5 * _indices[k] + 0.25 * (_indices[k - 1] + _indices[k + 1]);
                while (interval < last - 1 && index > interval + 1) interval++;
                double start = _nodes.GetFraction(interval);
                blocks = AddLevel(blocks, k, start + (index - interval) * (_nodes.GetFraction(interval + 1) - start));
            }
            SetLevels(blocks);
        }
        /// <summary>
        /// Where a fraction of the length lies among the nodes the edge was meshed with, node k being k.
        /// </summary>
        private double GetIndex(double fraction)
        {
            int last = _fractions.Length - 1;
            int interval = 0;
            while (interval < last - 1 && fraction > _nodes.GetFraction(interval + 1)) interval++;
            double start = _nodes.GetFraction(interval);
            return interval + (fraction - start) / (_nodes.GetFraction(interval + 1) - start);
        }
        /// <summary>
        /// Adds the new fraction of node k to the levels, pooling adjacent levels that are out of order, which projects
        /// the fractions onto the minimum spacing. Returns the number of levels.
        /// </summary>
        private int AddLevel(int blocks, int k, double fraction)
        {
            _levels[blocks] = fraction - k * _minimumInterval;
            _blockSizes[blocks] = 1;
            blocks++;
            while (blocks > 1 && _levels[blocks - 2] > _levels[blocks - 1])
            {
                blocks--;
                int size = _blockSizes[blocks - 1] + _blockSizes[blocks];
                _levels[blocks - 1] = (_levels[blocks - 1] * _blockSizes[blocks - 1] +
                                       _levels[blocks] * _blockSizes[blocks]) / size;
                _blockSizes[blocks - 1] = size;
            }
            return blocks;
        }
        private double SetLevels(int blocks)
        {
            int last = _fractions.Length - 1;
            double maxChange = 0;
            double upper = 1 - last * _minimumInterval;
            int node = 1;
            for (int b = 0; b < blocks; b++)
            {
                double level = Math.Max(0, Math.Min(upper, _levels[b]));
                for (int n = 0; n < _blockSizes[b]; n++, node++)
                {
                    double fraction = level + node * _minimumInterval;
                    maxChange = Math.Max(maxChange, Math.Abs(fraction - _fractions[node]));
                    _fractions[node] = fraction;
                }
            }
            return maxChange;
        }
        public void PrepareNodes(double tolerance)
        {
            for (int k = 1; k < _fractions.Length - 1; k++)
            {
                if (Math.Abs(_fractions[k] - _nodes.GetFraction(k)) <= tolerance)
                {
                    _fractions[k] = _nodes.GetFraction(k);
                    continue;
                }
                if (_proposedPoints == null)
                {
                    _proposedPoints = new double[_fractions.Length][];
                    _proposedParameters = new double[_fractions.Length];
                }
                _proposedPoints[k] = _nodes.EvaluateShellPoint(_fractions[k], out _proposedParameters[k]);
            }
        }
        public void ApplyNodes()
        {
            if (_proposedPoints == null) return;
            for (int k = 1; k < _fractions.Length - 1; k++)
            {
                if (_proposedPoints[k] != null) _nodes.SetShellPoint(k, _proposedPoints[k], _proposedParameters[k]);
            }
        }
        public void RestoreNodes()
        {
            if (_proposedPoints == null) return;
            for (int k = 1; k < _fractions.Length - 1; k++)
            {
                if (_proposedPoints[k] != null) _nodes.RestoreShellPoint(k);
            }
        }
        public void ResetNodes()
        {
            for (int k = 1; k < _fractions.Length - 1; k++) _fractions[k] = _nodes.GetFraction(k);
            _proposedPoints = null;
            _proposedParameters = null;
        }
    }
}
