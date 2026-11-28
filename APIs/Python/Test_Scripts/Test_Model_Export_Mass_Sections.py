#PrePoMax Python Script
#
#Checks the mass-section set names in successive live keyword previews and solver exports.
#
#  PrePoMax.com -f cube.inp -u MM_TON_S_C --script Test_Model_Export_Mass_Sections.py -g No -x Yes -w <dir>
#
#Expected process exit code 0.

import os


def preview():
    return "\n".join(entry.data for entry in pmx.model.keywords.get_model_keywords())


def check_mass_loads(deck, section):
    for load_type in ["Grav", "CENTRIF"]:
        assert section + "_All, " + load_type + "," in deck, section + " missing from " + load_type


surface = "MassSurface"
point = "PointMass"
distributed = "DistributedMass"
pmx.model.surfaces.add(surface, pmx.create.selection_by_ids('femodel', 'surfaces', 'surfaces', {1: [2]}))
pmx.model.sections.add_point_mass(point, pmx.create.selection_by_ids('femodel', 'nodes', 'nodes', [1]), 0.5)
pmx.model.sections.add_distributed_mass(distributed, pmx.create.region('SurfaceName', surface), 2)
pmx.model.steps.add_static("MassLoads")
loads = pmx.model.steps["MassLoads"].loads
for section in [point, distributed]:
    region = pmx.create.region('MassSection', section)
    loads.add_gravity("G" + section, region, f3=-9.81)
    loads.add_centrifugal("C" + section, region, axis_point=[0, 0, 0], axis_direction=[0, 0, 1], rotational_speed=50)

first = preview()
assert first == preview(), "Consecutive keyword previews differ"
for section in [point, distributed]:
    check_mass_loads(first, section)
for solver, export in [("Calculix", pmx.export_to_calculix), ("Abaqus", pmx.export_to_abaqus)]:
    path = os.path.join(pmx.work_directory, solver + "MassSections.inp")
    export(path)
    with open(path) as reader:
        deck = reader.read()
    for section in [point, distributed]:
        check_mass_loads(deck, section)

# An invalid section can still appear in the keyword editor. Its old generated set no longer exists.
pmx.model.surfaces.remove([surface])
updated = preview()
assert distributed + "_All" not in updated, "Keyword preview reused the distributed mass set from the preceding build"
check_mass_loads(updated, point)
print("TEST PASSED")
