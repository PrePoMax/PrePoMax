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
    /// The surface a mesh is smoothed and cleaned on in space (QuadSmoothing.SmoothRegionInSpace, QuadCleanup).
    /// </summary>
    public interface IQuadSurface
    {
        /// <summary>
        /// Moves the given nodes onto the surface and sets the normals of the surface there. A node the surface cannot
        /// take goes back to its previous position, given in the order of the nodes, and keeps its normal. Without
        /// previous positions the nodes are new ones, which stay where they are when the surface cannot take them.
        /// </summary>
        void Project(double[][] xyz, double[][] normals, IList<int> nodes, double[][] previousXyz);
    }
}
