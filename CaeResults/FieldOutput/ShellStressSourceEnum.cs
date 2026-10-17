// PrePoMax - Copyright (C) 2016-2026 Matej Borovinšek
//
// Licensed under the terms defined in the LICENSE file located in the root directory of this source code.
//
// Source code: https://gitlab.com/MatejB/PrePoMax
//
// Author: Matej Borovinšek
// Contributors:

using System;
using DynamicTypeDescriptor;

namespace CaeResults
{
    // The in-plane stresses of a ply: the stresses of the solver in the ply axes, or the plane stresses the elastic
    // constants of the ply give its in-plane strains, as classical lamination theory computes them
    [Serializable]
    public enum ShellStressSourceEnum
    {
        [StandardValue("Stress", DisplayName = "Stress field")]
        Stress,
        [StandardValue("Strain", DisplayName = "Strain field (plane stress)")]
        Strain
    }
}
