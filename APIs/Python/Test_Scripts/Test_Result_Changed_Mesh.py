#PrePoMax Python Script
#The PrePoMax object is available as 'pmx'
#
#Checks that the results of a mesh the model no longer has are not given the model parts. The model and its results
#share the hash name of the model, which a new mesh keeps, so the results are matched to the model mesh element by
#element. A fine mesh is solved and the part is meshed again coarser: the elements of the new mesh have the numbers of
#elements of the old one, but not their place, so the old results keep the one part of their own mesh instead of a
#model part rebuilt from the old elements that happen to have the new numbers next to a part of the rest, and the node
#sets of the model are not copied onto them. A mesh changed only to the other element order keeps its elements in
#place, so its results are still given the model part and its node sets.
#
#  PrePoMax.com -f L_bracket_r.step -u MM_TON_S_C --script Test_Result_Changed_Mesh.py -g No -x Yes -w <dir>
#
#Expected process exit code 0.
import os
import tempfile

PREFIX = "PyTest_Changed_Mesh_"


def check(condition, message):
    if not condition:
        raise Exception("TEST FAILED: " + message)


def mesh_part(part_name, max_size):
    geometry_part = pmx.model.geometry.parts[part_name]
    selection = pmx.create.selection_by_ids('geometry', 'geometry', 'parts', [geometry_part.id])
    setup_items = pmx.model.mesh.setup_items
    if (PREFIX + "Size") in list(setup_items.get_names()):
        setup_items.replace_meshing_parameters(PREFIX + "Size", max_size=max_size, min_size=max_size / 10.0)
    else:
        setup_items.add_meshing_parameters(PREFIX + "Size", selection, max_size=max_size, min_size=max_size / 10.0)
    report = pmx.model.mesh.generate([part_name])
    check(report[part_name] == "", "the part " + part_name + " was not meshed: " + report[part_name])


def solve(name, work):
    pmx.jobs.add(name, work_directory=work)
    status = pmx.jobs[name].run(timeout=300)
    check(status == "OK", "the analysis " + name + " did not succeed, got " + status + "\n" +
          pmx.jobs[name].get_output()[-3000:])
    return os.path.join(work, name + ".frd")


def open_results(file_name):
    pmx.results.clear()
    pmx.results.open(file_name)
    return pmx.results.current


def result_parts(result):
    parts = result.mesh.parts
    return {str(name): len(list(parts[name].get_element_ids())) for name in parts.get_names()}


def build_model(part_name):
    part_selection = pmx.create.selection_by_ids('femodel', 'elements', 'parts', [pmx.model.parts[part_name].id])
    pmx.model.materials.add(PREFIX + "Steel", pmx.create.elastic(210000, 0.3))
    pmx.model.sections.add_solid(PREFIX + "Section", PREFIX + "Steel", part_selection)
    node_ids = sorted(pmx.model.parts[part_name].get_node_ids())
    pmx.model.node_sets.add(PREFIX + "Fixed", pmx.create.selection_by_ids('femodel', 'nodes', 'nodes',
                                                                         node_ids[:len(node_ids) // 2]))
    step = pmx.model.steps.add_static(PREFIX + "Step")
    step.bcs.add_fixed(PREFIX + "Fix", pmx.create.region('NodeSetName', PREFIX + "Fixed"))


def test_remeshed_model(part_name, work):
    fine_file_name = solve(PREFIX + "Fine", work)
    fine_element_count = pmx.model.elements.count
    mesh_part(part_name, 20.0)
    coarse_element_count = pmx.model.elements.count
    check(coarse_element_count < fine_element_count,
          "the coarse mesh has no fewer elements than the fine one: " + str(coarse_element_count) + " of " +
          str(fine_element_count))
    result = open_results(fine_file_name)
    parts = result_parts(result)
    check(len(parts) == 1, "the results of the fine mesh were split into parts by the coarse model mesh: " + str(parts))
    check(list(parts.values())[0] == fine_element_count,
          "the part of the fine results does not hold all its elements: " + str(parts))
    node_set_names = [str(name).upper() for name in result.mesh.node_sets.get_names()]
    check((PREFIX + "Fixed").upper() not in node_set_names,
          "the node set of the coarse model mesh was copied onto the results of the fine mesh: " + str(node_set_names))


def test_changed_element_order(part_name, work):
    pmx.model.mesh.change_element_order([part_name], second_order=False)
    pmx.model.node_sets.remove([PREFIX + "Fixed"])
    node_ids = sorted(pmx.model.parts[part_name].get_node_ids())
    pmx.model.node_sets.add(PREFIX + "Fixed", pmx.create.selection_by_ids('femodel', 'nodes', 'nodes',
                                                                         node_ids[:len(node_ids) // 2]))
    linear_file_name = solve(PREFIX + "Linear", work)
    pmx.model.mesh.change_element_order([part_name], second_order=True)
    result = open_results(linear_file_name)
    parts = result_parts(result)
    check(list(parts.keys()) == [part_name],
          "the results of the linear mesh were not given the part of the parabolic model mesh: " + str(parts))
    node_set_names = [str(name).upper() for name in result.mesh.node_sets.get_names()]
    check((PREFIX + "Fixed").upper() in node_set_names,
          "the node set of the model was not copied onto the results of the same mesh of the other element order: " +
          str(node_set_names))


def main():
    work = tempfile.mkdtemp(prefix=PREFIX)
    part_names = list(pmx.model.geometry.parts.get_names())
    check(len(part_names) == 1, "the geometry is not one part: " + str(part_names))
    part_name = part_names[0]
    mesh_part(part_name, 5.0)
    pmx.model.mesh.change_element_order([part_name], second_order=False)
    build_model(part_name)
    test_remeshed_model(part_name, work)
    test_changed_element_order(part_name, work)
    print("TEST PASSED")


main()
