// PrePoMax - Copyright (C) 2016-2026 Matej Borovinšek
//
// Licensed under the terms defined in the LICENSE file located in the root directory of this source code.
//
// Source code: https://gitlab.com/MatejB/PrePoMax
//
// Author: Matej Borovinšek
// Contributors:

using CaeResults;

namespace APIs.Python
{
    public class ResultFieldOutputShellPlaneStressDraftAPI : ResultFieldOutputDraftAPI
    {
        // Variables
        private readonly string _fieldName;
        private readonly ShellStressSourceEnum _stressSource;


        // Properties
        [PythonAPI("The source tensor field name.", Returns = "str")]
        public string field_name { get { return _fieldName; } }
        [PythonAPI("Stress for the stresses of a stress field in the ply axes, or Strain for the plane stresses the " +
                   "elastic constants of the ply material give the in-plane strains of a strain field.", Returns = "str")]
        public string stress_source { get { return _stressSource.ToString(); } }


        // Constructor
        internal ResultFieldOutputShellPlaneStressDraftAPI(string name, string fieldName, ShellStressSourceEnum stressSource)
            : base(name)
        {
            _fieldName = fieldName;
            _stressSource = stressSource;
        }


        // Methods
        protected override ResultFieldOutput CreateFieldOutput()
        {
            ResultFieldOutputShellPlaneStress fieldOutput = new ResultFieldOutputShellPlaneStress(name, _fieldName);
            fieldOutput.StressSource = _stressSource;
            return fieldOutput;
        }
        public override string ToString()
        {
            return "FieldOutputShellPlaneStressDraft(name=" + StringListAPI.ToDisplayString(name) +
                   ", field_name=" + StringListAPI.ToDisplayString(_fieldName) +
                   ", stress_source=" + StringListAPI.ToDisplayString(_stressSource.ToString()) + ")";
        }
    }
}
