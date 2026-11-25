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
using System.Text;
using System.Threading.Tasks;
using PrePoMax;
using CaeModel;
using CaeMesh;
using CaeGlobals;
using PrePoMax.Python.Journal;
using Calculix = FileInOut.Output.Calculix;


namespace PrePoMax.Commands
{
    [PmxType(10549)]
    class CSetMaterials : PreprocessCommand
    {
        // Variables
        [PmxKey(0)] private string[] _oldMaterialNames;
        [PmxKey(1)] private Material[] _materials;
        [PmxKey(2)] private OrderedDictionary<int[], Calculix.CalculixUserKeyword> _userKeywords;
        [PmxKey(3)] private List<Calculix.CalculixKeywordOverride> _keywordOverrides;


        // Constructor
        public CSetMaterials(string[] oldMaterialNames, Material[] materials,
                             OrderedDictionary<int[], Calculix.CalculixUserKeyword> userKeywords,
                             List<Calculix.CalculixKeywordOverride> keywordOverrides)
            : base("Set materials")
        {
            _oldMaterialNames = oldMaterialNames;
            _materials = materials.DeepClone();
            _userKeywords = userKeywords?.DeepClone();
            _keywordOverrides = keywordOverrides?.DeepClone();
        }


        // Methods
        public override bool Execute(Controller receiver)
        {
            receiver.SetMaterials(_materials.DeepClone(), _userKeywords?.DeepClone(), _keywordOverrides?.DeepClone());
            return true;
        }
        public override void ToPython(PythonJournalWriter writer)
        {
            if (_oldMaterialNames.Length > 0)
                PythonModelItemEmitter.Remove(writer, "pmx.model.materials", _oldMaterialNames);
            foreach (var material in _materials) PythonMaterialEmitter.Add(writer, material);
            if (_userKeywords != null && _keywordOverrides != null)
                new CSetCalculixKeywordEditorState(_userKeywords, _keywordOverrides).ToPython(writer);
        }
        public override string GetCommandString()
        {
            return base.GetCommandString() + GetArrayAsString(_materials.Select(material => material.Name).ToArray());
        }
    }
}
