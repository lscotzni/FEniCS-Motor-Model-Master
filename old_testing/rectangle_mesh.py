import gmsh
import meshio
from math import sin, cos, tan, pi
import os

l = 2
w = 6

gmsh.initialize()
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("Rectangle") # Creating new model called motor
rect = gmsh.model.occ

lc = 5e-2 # target mesh size

rect.addPoint(0,0,0,lc,1)
rect.addPoint(w,0,0,lc,2)
rect.addPoint(w,l,0,lc,3)
rect.addPoint(0,l,0,lc,4)

''' Add next points for a second mesh to test different material properties '''
rect.addPoint(w,2*l,0,lc,5)
rect.addPoint(0,2*l,0,lc,6)


# rect.addLine(1,2,1)
# rect.addLine(2,4,2)
# rect.addLine(4,3,3)
# rect.addLine(3,1,4)
# ''' Add next points for a second mesh to test different material properties '''
# rect.addLine(2,5,5)
# rect.addLine(5,6,6)
# rect.addLine(6,4,7)

rect.addLine(1,2,1)
rect.addLine(2,3,2)
rect.addLine(3,4,3)
rect.addLine(4,1,4)
rect.addLine(3,5,5)
rect.addLine(5,6,6)
rect.addLine(6,4,7)

rect.addCurveLoop([1,2,3,4],1)
rect.addCurveLoop([5,6,7,-3],2)

rect.addPlaneSurface([1],1)
rect.addPlaneSurface([2],2)

# Assigning Groups
# gmsh.model.addPhysicalGroup(2,[1,2],100)
# gmsh.model.setPhysicalName(2,1,"Entire Domain")

rect.synchronize()

gmsh.model.addPhysicalGroup(2,[1],1)
gmsh.model.setPhysicalName(2,1,"Bottom")

gmsh.model.addPhysicalGroup(2,[2],2)
gmsh.model.setPhysicalName(2,2,"Top")

gmsh.model.addPhysicalGroup(1,[4],1)
gmsh.model.setPhysicalName(1,1,"B_Line")

gmsh.model.addPhysicalGroup(1,[6],2)
gmsh.model.setPhysicalName(1,2,"T_Line")

gmsh.model.addPhysicalGroup(1,[1,5],3)
gmsh.model.setPhysicalName(1,3,"L_Line")

gmsh.model.addPhysicalGroup(1,[3,7],4)
gmsh.model.setPhysicalName(1,4,"R_Line")

gmsh.model.addPhysicalGroup(1,[2],5)
gmsh.model.setPhysicalName(1,5,"M_Line")

gmsh.model.mesh.generate(dim=2)

gmsh.write("rect_mesh.msh")
gmsh.write("rect_mesh.vtk")

gmsh.finalize()

# ---------------------------------------------------

# filename = 'rect_mesh.vtk'

# mesh = meshio.read(
#     filename,
#     file_format = 'vtk'
# )
# points = mesh.points
# cells = mesh.cells
# meshio.write_points_cells(
#     'rect_mesh.xml',
#     points,
#     cells,
# )

# ---------------------------------------------------
# os.system('python3 msh2xdmf.py -d 2 2d_motor.msh')
# The msh2xdmf call to terminal window doesn't work here; errors 
# say that meshio and numpy are not available packages