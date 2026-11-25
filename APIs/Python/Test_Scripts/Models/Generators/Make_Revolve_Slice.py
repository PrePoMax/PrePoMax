#PrePoMax Python Script
# Writes revolve_slice.brep for Test_Model_Meshing_Revolve.py through GmshCommon (OCC kernel): a slice of a ball, the
# half disc of radius 33.5 in the plane y = 0 turned through 60 degrees about its diameter on the z axis. Every vertex
# of the part lies on that axis, so no edge leaves a vertex of its source along the revolution.
#
#  PrePoMax.com -u MM_TON_S_C --script Make_Revolve_Slice.py -g No -x Yes -w <work folder>
#  The file goes to the folder in PMX_MODEL_DIR, or to the Models folder above this script when it is not set.
import os
import sys
SRC = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "Make_Cap_Faces.py")
source = open(SRC).read()
exec(source[:source.index("Gmsh.Initialize()")])

RADIUS = 33.5

Gmsh.Initialize()
Gmsh.Model.Add("revolve_slice")
diameter = line(0.0, -RADIUS, 0.0, RADIUS)
rim = arc(0.0, 0.0, RADIUS, 90.0, -90.0)
half_disc = OCC.AddPlaneSurface(Array[Int32]([OCC.AddCurveLoop(Array[Int32]([diameter, rim]))]))
OCC.Revolve(dimtags(2, [half_disc]), 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, math.radians(60.0), None, no_ints(),
            no_doubles(), False)
OCC.Synchronize()
Gmsh.Write(os.path.join(OUT, "revolve_slice.brep"))
Gmsh.FinalizeAll()
print("Written revolve_slice.brep")
