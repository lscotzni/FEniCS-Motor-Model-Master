import gmsh
from math import sin, cos, tan, pi, floor, ceil
import numpy as np
import sys
from global_var import Ri, Ro_Ri, CCC_Ri, CCR_Ri

# Ri = 1              # Inner Radius
# Ro = 10 * Ri        # Outer Radius
# CCC = 6 * Ri        # Center of Current Circle
# CCR = 1 * Ri        # Radius of Current Circle

# Ri = 1                  # Inner Radius
Ro = Ro_Ri * Ri         # Outer Radius
CCC = CCC_Ri * Ri       # Center of Current Circle
CCR = CCR_Ri * Ri       # Radius of Current Circle

ks = 0.15           # Target Mesh Size

gmsh.initialize()
gmsh.option.setNumber("General.Terminal",1)

gmsh.model.add("Disk")
disk = gmsh.model.occ

# ----------------- ADDING POINTS -----------------
num_points = 0
# ORIGIN
disk.addPoint(0.0,0.0,0.0,ks,1)
num_points += 1
# INNER CIRCLE
disk.addPoint(Ri,0.0,0.0,ks,num_points+1)
num_points += 1
inner_points = []
inner_points.append([(0,num_points)])
for i in range(3):
    num_points += 1
    temp_in = disk.copy(inner_points[0])
    disk.rotate(temp_in,0,0,0,0,0,1,pi*(i+1)/2)
    inner_points.append(temp_in)

# OUTER CIRCLE
disk.addPoint(Ro,0.0,0.0,ks,num_points + 1)
num_points += 1
outer_points = []
outer_points.append([(0,num_points)])
for i in range(3):
    num_points += 1
    temp_out = disk.copy(outer_points[0])
    disk.rotate(temp_out,0,0,0,0,0,1,pi*(i+1)/2)
    outer_points.append(temp_out)



# UPPER CURRENT CIRCLE
disk.addPoint(0.0,CCC,0.0,ks,num_points + 1)
disk.addPoint(CCR,CCC,0.0,ks,num_points + 2)
num_points += 2

upper_circle_points = []
upper_circle_points.append([(0,num_points-1)])
upper_circle_points.append([(0,num_points)])
for i in range(3):
    temp = disk.copy(upper_circle_points[1])
    num_points += 1
    disk.rotate(temp,0,CCC,0,0,0,1,pi*(i+1)/2)
    upper_circle_points.append(temp)


# LOWER CURRENT CIRCLE
disk.addPoint(0.0,-CCC,0.0,ks,num_points + 1)
disk.addPoint(CCR,-CCC,0.0,ks,num_points + 2)
num_points += 2

lower_circle_points = []
lower_circle_points.append([(0,num_points-1)])
lower_circle_points.append([(0,num_points)])

for i in range(3):
    temp = disk.copy(lower_circle_points[1])
    num_points += 1
    disk.rotate(temp,0,-CCC,0,0,0,1,pi*(i+1)/2)
    lower_circle_points.append(temp)

print(inner_points[0])
print(outer_points)
print(upper_circle_points)
print(lower_circle_points)

# ----------------- CREATING LINES -----------------
# LINES DEFINING X-AXIS
num_lines = 0
disk.addLine(inner_points[0][0][1],outer_points[0][0][1],1)
disk.addLine(inner_points[2][0][1],outer_points[2][0][1],2)
num_lines += 2

# INNER CIRCLE
inner_circle_lines = []
for i in range(4):
    if i == len(inner_points) - 1:
        temp = disk.addCircleArc(inner_points[i][0][1],1,inner_points[0][0][1],num_lines + 1)
    else:
        temp = disk.addCircleArc(inner_points[i][0][1],1,inner_points[i+1][0][1],num_lines + 1)
    inner_circle_lines.append(temp)
    num_lines += 1

# OUTER CIRCLE
outer_circle_lines = []
for i in range(4):
    if i == len(outer_points) - 1:
        temp = disk.addCircleArc(outer_points[i][0][1],1,outer_points[0][0][1],num_lines + 1)
    else:
        temp = disk.addCircleArc(outer_points[i][0][1],1,outer_points[i+1][0][1],num_lines + 1)
    outer_circle_lines.append(temp)
    num_lines += 1

# UPPER CIRCLE
upper_circle_lines = []
for i in range(1,5):
    if i == len(upper_circle_points) - 1:
        temp = disk.addCircleArc(upper_circle_points[i][0][1],upper_circle_points[0][0][1],upper_circle_points[1][0][1],num_lines + 1)
    else:
        temp = disk.addCircleArc(upper_circle_points[i][0][1],upper_circle_points[0][0][1],upper_circle_points[i+1][0][1],num_lines + 1)
    upper_circle_lines.append(temp)
    num_lines += 1

# LOWER CIRCLE
lower_circle_lines = []
for i in range(1,5):
    if i == len(lower_circle_points) - 1:
        temp = disk.addCircleArc(lower_circle_points[i][0][1],lower_circle_points[0][0][1],lower_circle_points[1][0][1],num_lines + 1)
    else:
        temp = disk.addCircleArc(lower_circle_points[i][0][1],lower_circle_points[0][0][1],lower_circle_points[i+1][0][1],num_lines + 1)
    lower_circle_lines.append(temp)
    num_lines += 1

# ----------------- CREATING CURVE LOOPS -----------------
# UPPER HEMISPHERE CURVE LOOP
disk.addCurveLoop([1,outer_circle_lines[0],outer_circle_lines[1],2,-inner_circle_lines[1],-inner_circle_lines[0]],1)

# LOWER HEMISPHERE CURVE LOOP
disk.addCurveLoop([1,-outer_circle_lines[3],-outer_circle_lines[2],2,inner_circle_lines[2],inner_circle_lines[3]],2)

# UPPER CIRCLE CURVE LOOP
disk.addCurveLoop(upper_circle_lines,3)

# LOWER CIRCLE CURVE LOOP
disk.addCurveLoop(lower_circle_lines,4)

# ----------------- CREATING PLANE SURFACES -----------------
# UPPER HEMISPHERE PLANE SURFACE
disk.addPlaneSurface([1,3],1)

# LOWER HEMISPHERE PLANE SURFACE
disk.addPlaneSurface([2,4],2)

# UPPER CIRCLE PLANE SURFACE
disk.addPlaneSurface([3],3)

# LOWER CIRCLE PLANE SURFACE
disk.addPlaneSurface([4],4)

disk.synchronize()

# ----------------- SETTING PHYSICAL GROUPS -----------------
gmsh.model.addPhysicalGroup(2,[1],1)
gmsh.model.setPhysicalName(2,1,"Upper Hemisphere")

gmsh.model.addPhysicalGroup(2,[2],2)
gmsh.model.setPhysicalName(2,2,"Lower Hemisphere")

gmsh.model.addPhysicalGroup(2,[3],3)
gmsh.model.setPhysicalName(2,3,"Upper Wire")

gmsh.model.addPhysicalGroup(2,[4],4)
gmsh.model.setPhysicalName(2,4,"Lower Wire")

gmsh.model.addPhysicalGroup(1,inner_circle_lines,1)
gmsh.model.setPhysicalName(1,1,"Inner Boundary")

gmsh.model.addPhysicalGroup(1,outer_circle_lines,2)
gmsh.model.setPhysicalName(1,2,"Outer Boundary")

gmsh.model.mesh.generate(2)
gmsh.write("disk_mesh.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

import os
os.system('python3 msh2xdmf.py -d 2 disk_mesh.msh')