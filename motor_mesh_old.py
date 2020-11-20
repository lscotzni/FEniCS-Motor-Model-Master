import gmsh
from math import sin, cos, tan, pi, floor, ceil
import numpy as np
import os

'''
Mesh Generation of 2D Radial Flux PMSM; only half is modeled (180 degrees) due to periodicity
This version has too many subdomain definitions for FEniCS, but still creates a visible mesh
'''

p = 8 # pole pairs
s = 9 # number of stator slots for windings per 180 degrees
m = 3 # number of phases for stator winding current

gmsh.initialize() # gmsh must be initialized in Python before using functions
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("Motor_Old") # Creating new model called motor
motor = gmsh.model.occ

# gmsh.model.mesh.setOrder(2)

lc = 1e-2 # target mesh size

# Rotor Angles
theta_p = 2*pi/2/p
theta_m = .78 * theta_p
theta_b = .0595 * theta_p 
theta_g = (theta_p - theta_m - theta_b)/2 # angle of each piece on each side of magnet

# Stator Angles
theta_t = pi/s
theta_sso = .5 * theta_t
theta_ssi = .3 * theta_sso

Rr = .80
Rtm = .79
Rtb = .775
Rbb = .75
Rbm = .745
Rin = .60 # Inner Radius of Rotor

Rout = 1.15
Rsy = 1.03
Rssi = .83
Rs = .81

D = 2 # Domain Radius
Ras = (Rr+Rs)/2 # Radius splitting air-gap mesh of rotor and stator
RR = (Rin+Rr)/2 # Midpoint to cut air-gap mesh in Rotor
RS = (Rsy+Rout)/2 # # Midpoint to cut air-gap mesh in Stator


'''Rotor'''
# # Points
# Rotor Core
motor.addPoint(0,0,0,lc,1)
motor.addPoint(Rin,0,0,lc,2)
motor.addPoint(Rtm,0,0,lc,3)
motor.addPoint(Rtm*cos(theta_b/2),Rtm*sin(theta_b/2),0,lc,4)
motor.addPoint(Rr*cos(theta_b/2+theta_g),Rr*sin(theta_b/2+theta_g),0,lc,5)
motor.addPoint(Rr*cos(theta_b/2+theta_g+theta_m),Rr*sin(theta_b/2+theta_g+theta_m),0,lc,6)
motor.addPoint(Rtm*cos(theta_p-theta_b/2),Rtm*sin(theta_p-theta_b/2),0,lc,7)
motor.addPoint(Rtm*cos(theta_p),Rtm*sin(theta_p),0,lc,8)
motor.addPoint(Rin*cos(theta_p),Rin*sin(theta_p),0,lc,9)

# Magnet Slots
motor.addPoint(Rbb*cos(theta_b/2),Rbb*sin(theta_b/2),0,lc,10)
motor.addPoint(Rtb*cos(theta_b/2),Rtb*sin(theta_b/2),0,lc,11)
motor.addPoint(Rtm*cos(theta_b/2+theta_g),Rtm*sin(theta_b/2+theta_g),0,lc,12)
motor.addPoint(Rtm*cos(theta_b/2+theta_g+theta_m),Rtm*sin(theta_b/2+theta_g+theta_m),0,lc,13)
motor.addPoint(Rtb*cos(theta_p-theta_b/2),Rtb*sin(theta_p-theta_b/2),0,lc,14)
motor.addPoint(Rbb*cos(theta_p-theta_b/2),Rbb*sin(theta_p-theta_b/2),0,lc,15)
motor.addPoint(Rbb*cos(theta_b/2+theta_g+theta_m),Rbb*sin(theta_b/2+theta_g+theta_m),0,lc,16)
motor.addPoint(Rbm*cos(theta_b/2+theta_g+theta_m),Rbm*sin(theta_b/2+theta_g+theta_m),0,lc,17)
motor.addPoint(Rbm*cos(theta_b/2+theta_g),Rbm*sin(theta_b/2+theta_g),0,lc,18)
motor.addPoint(Rbb*cos(theta_b/2+theta_g),Rbb*sin(theta_b/2+theta_g),0,lc,19)

# # Connecting Lines
# Rotor Core
motor.addLine(2,3,1)
motor.addLine(3,4,2)
motor.addLine(4,5,3)
motor.addCircleArc(5,1,6,4)
motor.addLine(6,7,5)
motor.addLine(7,8,6)
motor.addLine(8,9,7)
motor.addCircleArc(9,1,2,8)

# Magnet Slots
motor.addLine(10,11,9)
motor.addLine(11,12,10)
motor.addCircleArc(12,1,13,11)
motor.addLine(13,14,12)
motor.addLine(14,15,13)
motor.addCircleArc(15,1,16,14)
motor.addLine(16,17,15)
motor.addCircleArc(17,1,18,16)
motor.addLine(18,19,17)
motor.addCircleArc(19,1,10,18)

# Magnets
motor.addLine(18,12,19)
''' Next is Line 11'''
motor.addLine(13,17,21)
''' Next is Line 16'''

# Outer Magnet Piece
''' Line 9 '''
''' Line 10 '''
motor.addLine(12,19,22)
''' Line 18 '''

''' Line 13 '''
''' Line 14 '''
motor.addLine(16,13,23)
''' Line 12 '''

# Curved Loop
motor.addCurveLoop([1,2,3,4,5,6,7,8],1)
motor.addCurveLoop([9,10,11,12,13,14,15,16,17,18],2)
motor.addCurveLoop([19,11,21,16],3)
motor.addCurveLoop([9,10,22,18],4)
motor.addCurveLoop([13,14,23,12],5)

# Plane Surface
motor.addPlaneSurface([1,2],1)
motor.addPlaneSurface([3],2)
motor.addPlaneSurface([4],3)
motor.addPlaneSurface([5],4)

'''Stator'''
# # Points
# Core
motor.addPoint(Rsy,0,0,lc,20)
motor.addPoint(Rout,0,0,lc,21)
motor.addPoint(Rout*cos(theta_t),Rout*sin(theta_t),0,lc,22)
motor.addPoint(Rsy*cos(theta_t),Rsy*sin(theta_t),0,lc,23)
motor.addPoint(Rsy*cos(theta_t-theta_sso/2),Rsy*sin(theta_t-theta_sso/2),0,lc,24)
motor.addPoint(Rssi*cos(theta_t-theta_sso/2),Rssi*sin(theta_t-theta_sso/2),0,lc,25)
motor.addPoint(Rssi*cos(theta_t-theta_ssi/2),Rssi*sin(theta_t-theta_ssi/2),0,lc,26)
motor.addPoint(Rs*cos(theta_t-theta_ssi/2),Rs*sin(theta_t-theta_ssi/2),0,lc,27)
motor.addPoint(Rs*cos(theta_ssi/2),Rs*sin(theta_ssi/2),0,lc,28)
motor.addPoint(Rssi*cos(theta_ssi/2),Rssi*sin(theta_ssi/2),0,lc,29)
motor.addPoint(Rssi*cos(theta_sso/2),Rssi*sin(theta_sso/2),0,lc,30)
motor.addPoint(Rsy*cos(theta_sso/2),Rsy*sin(theta_sso/2),0,lc,31)

# Windings
motor.addPoint(Rssi,0,0,lc,32)
motor.addPoint(Rssi*cos(theta_t),Rssi*sin(theta_t),0,lc,33)

# # Connecting Lines
# Core
motor.addLine(20,21,24)
motor.addCircleArc(21,1,22,25)
motor.addLine(22,23,26)
motor.addCircleArc(23,1,24,27)
motor.addLine(24,25,28)
motor.addCircleArc(25,1,26,29)
motor.addLine(26,27,30)
motor.addCircleArc(27,1,28,31)
motor.addLine(28,29,32)
motor.addCircleArc(29,1,30,33)
motor.addLine(30,31,34)
motor.addCircleArc(31,1,20,35)

# Windings
motor.addLine(32,20,36)
# Line -35
# Line -34
motor.addCircleArc(30,1,32,37)
motor.addLine(23,33,38)
motor.addCircleArc(33,1,25,39)
# Line -28
# Line -27

# Curved Loop
motor.addCurveLoop([24,25,26,27,28,29,30,31,32,33,34,35],6) # Core
motor.addCurveLoop([36,-35,-34,37],7) # Winding
motor.addCurveLoop([38,39,-28,-27],8) # Winding

# Plane Surface
motor.addPlaneSurface([6],5)
motor.addPlaneSurface([7],6)
motor.addPlaneSurface([8],7)


'''Air-gap'''
'''Rotor Air Gap'''
# Points
# Need point 1 for origin (Inner)
motor.addPoint(RR,0,0,lc,34) # Inner/Outer
motor.addPoint(RR*cos(theta_p),RR*sin(theta_p),0,lc,35) # Inner/Outer
motor.addPoint(Ras,0,0,lc,36) # Outer
motor.addPoint(Ras*cos(theta_p),Ras*sin(theta_p),0,lc,37) # Outer

# Connecting Lines
motor.addLine(1,34,40) # Inner
motor.addCircleArc(34,1,35,41) # Inner/Outer
motor.addLine(35,1,42) # Inner
motor.addLine(34,36,43)
motor.addCircleArc(36,1,37,44)
motor.addLine(37,35,45)

# Curve Loop
motor.addCurveLoop([40,41,42],9)
motor.addCurveLoop([43,44,45,-41],10)

# Plane Surface
motor.addPlaneSurface([9],8)
motor.addPlaneSurface([10],9)

'''Stator Air Gap'''
# Points
# Point 36 is needed for stator air gap
motor.addPoint(RS,0,0,lc,38)
motor.addPoint(RS*cos(theta_t),RS*sin(theta_t),0,lc,39)
motor.addPoint(Ras*cos(theta_t),Ras*sin(theta_t),0,lc,40)
motor.addPoint(D,0,0,lc,41)
motor.addPoint(D*cos(theta_t),D*sin(theta_t),0,lc,42)

# Connecting Lines
motor.addLine(36,38,46)
motor.addCircleArc(38,1,39,47)
motor.addLine(39,40,48)
motor.addCircleArc(40,1,36,49)

motor.addLine(38,41,50)
motor.addCircleArc(41,1,42,51)
motor.addLine(42,39,52)
# Line -47

# Curve Loop
motor.addCurveLoop([46,47,48,49],11)
motor.addCurveLoop([50,51,52,-47],12)

# Plane Surface
motor.addPlaneSurface([11],10)
asdf = motor.addPlaneSurface([12],11)

'''Whole Boundary'''
# Point 
motor.addPoint(-D,0,0,lc,100)
# Lines
motor.addLine(100,41,100)
motor.addCircleArc(41,1,100,101)
# Curve Loop
motor.addCurveLoop([100,101],100)



# ---------------------------- Boolean Operations For Air Gap ----------------------------
# '''Rotor'''
motor.cut([(2,8)],[(2,1)],12,removeObject=True,removeTool=False) # Inner Air Gap Cut
motor.cut([(2,9)],[(2,1),(2,2),(2,3),(2,4)],13,removeObject=True,removeTool=False) # Outer Air Gap Cut

# '''Stator'''
motor.cut([(2,10)],[(2,5),(2,6),(2,7)],14,removeObject=True,removeTool=False) # Inner Air Gap Cut
motor.cut([(2,11)],[(2,5)],15,removeObject=True,removeTool=False) # Outer Air Gap Cut

# Outline the solid portion in a color like red/black, airgaps with another

# --------------------------------------- Rotation ---------------------------------------
'''Rotor Rotation'''
rgeom = [0 for i in range(p-1)] # Rotor Geometry Copy
rcrot = [0 for i in range(p-1)] # Rotor Core Rotation
rmrot = [0 for i in range(p-1)] # Rotor Magnet Rotation
rrerot = [0 for i in range(p-1)] # Rotor Right Edge Rotation
rlerot = [0 for i in range(p-1)] # Rotor Left Edge Rotation
riarot = [0 for i in range(p-1)] # Rotor Inner Air Gap Rotation
roarot = [0 for i in range(p-1)] # Rotor Outer Air Gap Rotation

for n in range(0,p-1):
    rgeom[n] = motor.copy([(2,1),(2,2),(2,3),(2,4),(2,12),(2,13)])
    rcrot[n] = motor.rotate([(rgeom[n][0])],0,0,0,0,0,1,(n+1) * pi/p)
    rmrot[n] = motor.rotate([(rgeom[n][1])],0,0,0,0,0,1,(n+1) * pi/p)
    rrerot[n] = motor.rotate([(rgeom[n][2])],0,0,0,0,0,1,(n+1) * pi/p)
    rlerot[n] = motor.rotate([(rgeom[n][3])],0,0,0,0,0,1,(n+1) * pi/p)
    riarot[n] = motor.rotate([(rgeom[n][4])],0,0,0,0,0,1,(n+1) * pi/p)
    roarot[n] = motor.rotate([(rgeom[n][5])],0,0,0,0,0,1,(n+1) * pi/p)

'''Stator Rotation'''
sgeom = [0 for i in range(s-1)] # Stator Geometry Copy
scrot = [0 for i in range(s-1)] # Stator Core Rotation
srwrot = [0 for i in range(s-1)] # Stator Right Winding Rotation
slwrot = [0 for i in range(s-1)] # Stator Left Winding Rotation
siarot = [0 for i in range(s-1)] # Stator Inner Air Gap Rotation
soarot = [0 for i in range(s-1)] # Stator Outer Air Gap Rotation

for n in range(0,s-1):
    sgeom[n] = motor.copy([(2,5),(2,6),(2,7),(2,14),(2,15)])
    scrot[n] = motor.rotate([(sgeom[n][0])],0,0,0,0,0,1,(n+1) * pi/s)
    srwrot[n] = motor.rotate([(sgeom[n][1])],0,0,0,0,0,1,(n+1) * pi/s)
    slwrot[n] = motor.rotate([(sgeom[n][2])],0,0,0,0,0,1,(n+1) * pi/s)
    siarot[n] = motor.rotate([(sgeom[n][3])],0,0,0,0,0,1,(n+1) * pi/s)
    soarot[n] = motor.rotate([(sgeom[n][4])],0,0,0,0,0,1,(n+1) * pi/s)

''' Synchronize Mesh '''
motor.synchronize()
# -------------------------------------------------------------------------------------------

# ''' Physical Groups for FEniCS '''
# # Creating numpy array for Rotor and Stator Geometry
# rgeom_array = np.array(rgeom,dtype = int)
# sgeom_array = np.array(sgeom,dtype = int)

# '''Air Gap / Domain Definition'''
# air_gap_group = [0 for i in range(2*(2*p + s)-1)]
# air_gap_group[0:5] = [12,13,14,15,3,4]

# for n in range(0,p-1): # for Rotor 
#     air_gap_group[6+2*n] = rgeom_array[n,4,1]
#     air_gap_group[6+2*n+1] = rgeom_array[n,5,1]

# for n in range(0,s-1): # for Stator
#     air_gap_group[6+2*(p+n-1)] = sgeom_array[n,3,1]
#     air_gap_group[6+2*(p+n-1)+1] = sgeom_array[n,4,1]

# for n in range(0,p-1): # for Rotor Bridge Air Gaps
#     air_gap_group[6+2*(2*p+n-1)] = rgeom_array[n,2,1]
#     air_gap_group[6+2*(2*p+n-1)+1] = rgeom_array[n,3,1]

# gmsh.model.addPhysicalGroup(2,air_gap_group,1)
# gmsh.model.setPhysicalName(2, 1, "Air Gap")

# # -------------------------------------------------------------------------------------

# # '''TESTING BOUNDARY REQUIREMENTS FOR XDMF CONVERSION'''
# # air_gap_group[0] = 15
# # for n in range(1):
# #     air_gap_group[n+1] = sgeom_array[n,4,1]

# # gmsh.model.addPhysicalGroup(2,[air_gap_group[0]],1)
# # gmsh.model.setPhysicalName(2, 1, "Air Gap 1")

# # gmsh.model.addPhysicalGroup(2,[air_gap_group[1]],2)
# # gmsh.model.setPhysicalName(2, 2, "Air Gap 2")

# # gmsh.model.addPhysicalGroup(1,[25 ,-112 ,-113 ,-114],1000)
# # gmsh.model.setPhysicalName(1,1000,'Pizza slice 1')

# # gmsh.model.addPhysicalGroup(1,[749,-750,-751,-752],1001)
# # gmsh.model.setPhysicalName(1,1001,'Pizza slice 2')

# # aaa = gmsh.model.getBoundary([(2,15)],combined=False) # boundary of first slice of stator airgap 
# # bbb = gmsh.model.getBoundary([(2,62)],combined=False) # boundary of second slice of stator airgap

# # -------------------------------------------------------------------------------------

# '''Rotor Core'''
# rotor_core_group = [0 for i in range(p)]
# rotor_core_group[0] = 1

# for n in range(0,p-1):
#     rotor_core_group[1+n] = rgeom_array[n,0,1]

# gmsh.model.addPhysicalGroup(2,rotor_core_group,2)
# gmsh.model.setPhysicalName(2, 2, "Rotor Core")

# '''Stator Core'''
# stator_core_group = [0 for i in range(s)]
# stator_core_group[0] = 5

# for n in range(0,s-1):
#     stator_core_group[1+n] = sgeom_array[n,0,1]

# gmsh.model.addPhysicalGroup(2,stator_core_group,3)
# gmsh.model.setPhysicalName(2, 3, "Stator Core")

# '''Magnets'''
# # There are two poles here, so we will do alternating groups; discrepancy appears for odd number of poles

# # Odd Magnets (starting from positive x-axis)
# magnet_group_1 = [0 for i in range(ceil(p/2))]
# magnet_group_1[0] = 2

# for n in range(0,ceil(p/2)-1):
#     magnet_group_1[1+n] = rgeom_array[2*n+1,1,1]

# gmsh.model.addPhysicalGroup(2,magnet_group_1,4)
# gmsh.model.setPhysicalName(2, 4, "Odd Magnets")

# # Even Magnets (starting from positive x-axis)
# if (p%2 != 0):
#     magnet_group_2 = [0 for i in range(floor(p/2))]
#     for n in range(floor(p/2)):
#         magnet_group_2[n] = rgeom_array[2*n,1,1]
# else:
#     magnet_group_2 = [0 for i in range(ceil(p/2))]
#     for n in range(ceil(p/2)):
#         magnet_group_2[n] = rgeom_array[2*n,1,1]

# gmsh.model.addPhysicalGroup(2,magnet_group_2,5)
# gmsh.model.setPhysicalName(2, 5, "Even Magnets")

# '''Stator Windings'''
# # Label each phase as A, B, and C; group each into 1 (+) and 2 (-)

# stator_winding_A1 = [6,60,64]
# gmsh.model.addPhysicalGroup(2,stator_winding_A1,6)
# gmsh.model.setPhysicalName(2, 6, "Stator Phase A+")

# stator_winding_A2 = [7,59,65]
# gmsh.model.addPhysicalGroup(2,stator_winding_A2,7)
# gmsh.model.setPhysicalName(2, 7, "Stator Phase A-")

# stator_winding_B1 = [69,75,79]
# gmsh.model.addPhysicalGroup(2,stator_winding_B1,8)
# gmsh.model.setPhysicalName(2, 8, "Stator Phase B+")

# stator_winding_B2 = [70,74,80]
# gmsh.model.addPhysicalGroup(2,stator_winding_B2,9)
# gmsh.model.setPhysicalName(2, 9, "Stator Phase B-")

# stator_winding_C1 = [84,90,94]
# gmsh.model.addPhysicalGroup(2,stator_winding_C1,10)
# gmsh.model.setPhysicalName(2, 10, "Stator Phase C+")

# stator_winding_C2 = [85,89,95]
# gmsh.model.addPhysicalGroup(2,stator_winding_C2,11)
# gmsh.model.setPhysicalName(2, 11, "Stator Phase C-")

# s_w_C2_b = [0 for i in range(len(stator_winding_C2))]
# for n in range(len(s_w_C2_b)):
#     s_w_C2_b[n] = gmsh.model.getBoundary([(2,stator_winding_C2[n])])

# gmsh.model.addPhysicalGroup(1,s_w_C2_b,100)
# gmsh.model.setPhysicalName(1, 100, "Stator Phase C- Boundary")

# '''Physical Group for Boundary'''
# gmsh.model.addPhysicalGroup(1,[100, 101],1)
# gmsh.model.setPhysicalName(1,1,"Domain Boundary")

# Need to add boundary of domain for approaching zero and and x-axis for periodic BCs

# -------------------------------------------------------------------------------------------
# aaa = gmsh.model.getPhysicalGroups(1)


gmsh.model.mesh.generate(2)

# gmsh.model.mesh.recombine()

gmsh.write("old_motor_mesh.msh")
# gmsh.write("2d_motor_test.msh")
# gmsh.write("2d_motor.vtk")

# print(rgeom_array[:p, 4, 1])
# print(rgeom_array[:p, 5, 1])
# print(sgeom_array[:p, 3, 1])
# print(sgeom_array[:p, 4, 1])
# print(air_gap_group)
# print(asdf)
# print(aaa)


# print(aaa)
# print(air_gap_group)
# print(bbb)

gmsh.finalize()



# -------------------------------------------------------------------------------------------
# '''Converting to XML file for FEniCS'''

# filename = '2d_motor.vtk'

# mesh = meshio.read(
#     filename,
#     file_format = 'vtk'
# )
# points = mesh.points
# cells = mesh.cells
# meshio.write_points_cells(
#     '2d_motor.xml',
#     points,
#     cells,
# )

# -------------------------------------------------------------------------------------------

# msh = meshio.read("2d_motor.msh")

# # meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
# # meshio.write("mf.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
# #                                     cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))

# for cell in msh.cells:
#     if cell.type == "triangle":
#         triangle_cells = cell.data
#     elif  cell.type == "line":
#         line_cells = cell.data

# for key in msh.cell_data_dict["gmsh:physical"].keys():
#     if key == "line":
#         line_data = msh.cell_data_dict["gmsh:physical"][key]
#     elif key == "triangle":
#         triangle_data = msh.cell_data_dict["gmsh:physical"][key]

# triangle_mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells})
# line_mesh =meshio.Mesh(points=msh.points,
#                            cells=[("line", line_cells)],
#                            cell_data={"name_to_read":[line_data]})
# meshio.write("2d_motor.xdmf", triangle_mesh)

# meshio.write("2d_motor_mf.xdmf", line_mesh)

# -------------------------------------------------------------------------------------------
# meshio conversion this way creates a working xdmf file, but
# the mesh has mixed shapes which can't be read by FEniCS
# os.system('meshio-convert 2d_motor.msh 2d_motor_terminal.xdmf')


# os.system('python3 msh2xdmf.py -d 2 old_motor_mesh.msh')
# THIS LINE WORKS only when called using the general terminal window.
# The msh2xdmf call to terminal window in VSCode doesn't work here; errors 
# say that meshio and numpy are not available packages

'''
NOTES

using msh2xdmf via direct terminal on the computer works, but there is a strange error with
|gmsh:bounding_entities| = 2 in the table, which randomly showed up after trying to download 
meshio via conda in both python and bash terminal windows in VSCode. The solution might be to
delete these downloads from those terminals. 




'''