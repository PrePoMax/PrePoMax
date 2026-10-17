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
    public class ResultFieldOutputCompositeFailureDraftAPI : ResultFieldOutputDraftAPI
    {
        // Variables
        private readonly string _fieldName;
        private readonly CompositeFailureCriterionEnum _criterion;
        private readonly ShellStressSourceEnum _stressSource;


        // Properties
        [PythonAPI("The source tensor field name.", Returns = "str")]
        public string field_name { get { return _fieldName; } }
        [PythonAPI("The failure criterion: MaxStress, TsaiHill, TsaiWu or Hashin.", Returns = "str")]
        public string criterion { get { return _criterion.ToString(); } }
        [PythonAPI("Stress for the stresses of a stress field in the ply axes, or Strain for the plane stresses the " +
                   "elastic constants of the ply material give the in-plane strains of a strain field.", Returns = "str")]
        public string stress_source { get { return _stressSource.ToString(); } }


        // Constructor
        internal ResultFieldOutputCompositeFailureDraftAPI(string name, string fieldName, CompositeFailureCriterionEnum criterion,
                                                           ShellStressSourceEnum stressSource)
            : base(name)
        {
            _fieldName = fieldName;
            _criterion = criterion;
            _stressSource = stressSource;
        }


        // Methods
        protected override ResultFieldOutput CreateFieldOutput()
        {
            ResultFieldOutputCompositeFailure fieldOutput = new ResultFieldOutputCompositeFailure(name, _fieldName, _criterion);
            fieldOutput.StressSource = _stressSource;
            return fieldOutput;
        }
        public override string ToString()
        {
            return "FieldOutputCompositeFailureDraft(name=" + StringListAPI.ToDisplayString(name) +
                   ", field_name=" + StringListAPI.ToDisplayString(_fieldName) +
                   ", criterion=" + StringListAPI.ToDisplayString(_criterion.ToString()) +
                   ", stress_source=" + StringListAPI.ToDisplayString(_stressSource.ToString()) + ")";
        }
    }
}
