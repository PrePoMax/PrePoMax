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

namespace CaeMesh
{
    /// <summary>
    /// The mappings available for the source faces of one sweep. A submapped face may keep a lens mapping in
    /// reserve until its edge counts or its trial mesh decide which mapping to use.
    /// </summary>
    internal sealed class GmshSweepSources
    {
        // Variables                                                                                                                
        private readonly Dictionary<int, TransfiniteFace> _structured = new Dictionary<int, TransfiniteFace>();
        private readonly Dictionary<int, SubmapFace> _submaps = new Dictionary<int, SubmapFace>();
        private readonly Dictionary<int, RingFace> _rings = new Dictionary<int, RingFace>();
        private readonly Dictionary<int, LensFace> _lenses = new Dictionary<int, LensFace>();
        private readonly Dictionary<int, LensFace> _fallbackLenses = new Dictionary<int, LensFace>();
        private readonly HashSet<int> _thinLenses = new HashSet<int>();
        private readonly HashSet<int> _raisedRoundHoles = new HashSet<int>();
        private bool _edgeTrialsStarted;


        // Properties                                                                                                               
        public IReadOnlyDictionary<int, TransfiniteFace> Structured { get { return _structured; } }
        public IReadOnlyDictionary<int, SubmapFace> Submaps { get { return _submaps; } }
        public IReadOnlyDictionary<int, RingFace> Rings { get { return _rings; } }
        public IReadOnlyDictionary<int, LensFace> Lenses { get { return _lenses; } }
        public IReadOnlyDictionary<int, LensFace> FallbackLenses { get { return _fallbackLenses; } }


        // Methods                                                                                                                  
        public void Register(GmshSurface surface, TransfiniteFace structured, SubmapFace submap, RingFace ring,
                             LensFace lens, Func<RingFace, bool> isThinRing)
        {
            int id = surface.Id;
            surface.Transfinite = structured != null && structured.NumberOfSides == 4 && surface.FourSided &&
                                  structured.SideEdgeIds.All(side => side.Length == 1);
            if (structured != null)
            {
                if (!surface.Transfinite) _structured.Add(id, structured);
                return;
            }
            // A thin ring takes precedence over a submap. Only test the ring thickness when both exist.
            if (ring != null && (submap == null || isThinRing(ring)))
            {
                _rings.Add(id, ring);
                return;
            }
            if (submap != null && (lens == null || !GmshLensMesher.HasNoStraightEdge(lens)))
            {
                _submaps.Add(id, submap);
                if (lens != null) _fallbackLenses.Add(id, lens);
                return;
            }
            if (lens != null) _lenses.Add(id, lens);
        }
        public bool BeginEdgeTrials()
        {
            if (_edgeTrialsStarted) return false;
            _edgeTrialsStarted = true;
            return true;
        }
        public void MarkThinLenses(Dictionary<int, int> naturalCounts)
        {
            foreach (var face in _lenses.Values.Concat(_fallbackLenses.Values))
            {
                if (GmshLensMesher.IsThin(face, naturalCounts)) _thinLenses.Add(face.FaceId);
            }
        }
        public bool TryUseLens(int id)
        {
            LensFace lens;
            if (!_fallbackLenses.TryGetValue(id, out lens)) return false;
            _fallbackLenses.Remove(id);
            _submaps.Remove(id);
            _lenses.Add(id, lens);
            return true;
        }
        public bool Contains(int id)
        {
            return _structured.ContainsKey(id) || _submaps.ContainsKey(id) || _rings.ContainsKey(id) || _lenses.ContainsKey(id);
        }
        public bool IsThinLens(int id) { return _thinLenses.Contains(id); }
        public bool HasRaisedRoundHole(int id) { return _raisedRoundHoles.Contains(id); }
        public void MarkRaisedRoundHole(int id) { _raisedRoundHoles.Add(id); }
        public void RemoveStructured(int id) { _structured.Remove(id); }
        public void RemoveRing(int id) { _rings.Remove(id); }
        public void RemoveLens(int id) { _lenses.Remove(id); }
        public void RemoveSubmap(int id) { _submaps.Remove(id); }
    }
}
