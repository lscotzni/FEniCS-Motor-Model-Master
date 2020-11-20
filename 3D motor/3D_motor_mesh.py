import gmsh
from math import sin, cos, tan, pi
import matplotlib.pyplot as plt
''' Mesh generates geometry of YASA PMSM in the axial form for a given number of poles.
The rotor disks and stator area are generated separately
'''

# NOTE: do all extruding at the end; extruding process uses arbitrary points and the
# point tags will not work if extruded beforehand

p = 8 # pole pairs
t = 15 # number of stator teeth
R = 4. # outer radius
r = 2.5 # inner radius
rl = 1. # rotor disk thickness/length
sl = 3. # stator thickness/length
mt = .1 # magnet thickness
wt = .25 # stator winding thickness
ag = 1. # air gap between magnets and stator
D = 10. # domain size

y2 = .4 # depth of angled portion of stator core tooth tip
y3 = .2 # length of straight portion of stator core tooth tip (closest to magnets)

l = sl + 2 * (mt + ag + y2 + y3)

sm = mt + ag + y2 + y3 + sl/2 # middle of stator relative to origin

gmsh.initialize() # gmsh must be initialized in Python before using functions
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("Motor") # Creating new model called motor
motor = gmsh.model.geo

lc = 1e-1 # target mesh size

# Magnet Angular Sweep
angle_m = 2*pi/(2*p)/2
sweep_m = .8 * angle_m

# Stator Tooth Angular Sweep
angle_s = 2*pi/t
sweep_s = .8 * angle_s


'''Back of the Rotor'''
# Points
motor.addPoint(0,0,0,lc,1) # center of circle for back of rotor
motor.addPoint(R,0,0,lc,2) # Outer Radius
motor.addPoint(0,R,0,lc,3)
motor.addPoint(-R,0,0,lc,4)
motor.addPoint(0,-R,0,lc,5)
motor.addPoint(r,0,0,lc,6) # Inner Radius
motor.addPoint(0,r,0,lc,7)
motor.addPoint(-r,0,0,lc,8)
motor.addPoint(0,-r,0,lc,9)
# Circular Arcs
motor.addCircleArc(2,1,3,1) # Outer
motor.addCircleArc(3,1,4,2)
motor.addCircleArc(4,1,5,3)
motor.addCircleArc(5,1,2,4)
motor.addCircleArc(6,1,7,5) # Inner
motor.addCircleArc(7,1,8,6)
motor.addCircleArc(8,1,9,7)
motor.addCircleArc(9,1,6,8)
# Curved Loop
motor.addCurveLoop([1, 2, 3, 4], 1)
motor.addCurveLoop([5, 6, 7, 8], 2)
# Plane Surface
motor.addPlaneSurface([1, 2], 1)

'''Front of the Rotor'''
# Points
motor.addPoint(0,0,l,lc,10) # center of circle for front of rotor
motor.addPoint(R,0,l,lc,11) # Outer Radius
motor.addPoint(0,R,l,lc,12)
motor.addPoint(-R,0,l,lc,13)
motor.addPoint(0,-R,l,lc,14)
motor.addPoint(r,0,l,lc,15) # Inner Radius
motor.addPoint(0,r,l,lc,16)
motor.addPoint(-r,0,l,lc,17)
motor.addPoint(0,-r,l,lc,18)
# Circular Arcs
motor.addCircleArc(11,10,12,9)
motor.addCircleArc(12,10,13,10)
motor.addCircleArc(13,10,14,11)
motor.addCircleArc(14,10,11,12)
motor.addCircleArc(15,10,16,13)
motor.addCircleArc(16,10,17,14)
motor.addCircleArc(17,10,18,15)
motor.addCircleArc(18,10,15,16)
# Curved Loop
motor.addCurveLoop([9, 10, 11, 12], 3)
motor.addCurveLoop([13, 14, 15, 16], 4)
# Plane Surface
motor.addPlaneSurface([3, 4], 2)

'''Stator Winding'''
Roww = R * tan(sweep_s/2) # Outer Radius outside winding area width (from vertical centerline)
roww = r * tan(sweep_s/2) # Inner Radius outside winding area width (from vertical centerline)
Riww = (R - wt - wt/sin(sweep_s/2)) * tan(sweep_s/2) # Outer Radius inside winding area width (from vertical centerline)
riww = (r + wt - wt/sin(sweep_s/2)) * tan(sweep_s/2) # Inner Radius inside winding area width (from vertical centerline)

# Points
motor.addPoint(0,0,sm - sl/2,lc,19)
# motor.addPoint(0,R,sm,lc,) # can remove
motor.addPoint(Roww,R,sm - sl/2,lc,20)
motor.addPoint(-Roww,R,sm - sl/2,lc,21)
# motor.addPoint(0,r,sm,lc,) # can remove
motor.addPoint(roww,r,sm - sl/2,lc,22)
motor.addPoint(-roww,r,sm - sl/2,lc,23)
# motor.addPoint(0,R-wt,sm,lc,) # can remove
motor.addPoint(Riww,R-wt,sm - sl/2,lc,24)
motor.addPoint(-Riww,R-wt,sm - sl/2,lc,25)
# motor.addPoint(0,r+wt,sm - sl/2,lc,) # can remove
motor.addPoint(riww,r+wt,sm - sl/2,lc,26)
motor.addPoint(-riww,r+wt,sm - sl/2,lc,27)

# Connecting Lines
aaa = motor.addLine(20,21,17)
motor.addLine(21,23,18)
motor.addLine(23,22,19)
motor.addLine(22,20,20)

motor.addLine(24,25,21)
motor.addLine(25,27,22)
motor.addLine(27,26,23)
motor.addLine(26,24,24)

# Curved Loop
bbb = motor.addCurveLoop([17,18,19,20],5)
motor.addCurveLoop([21,22,23,24],6)

# Plane Surface
ccc = motor.addPlaneSurface([5,6],3)

# print(aaa)
# print(bbb)
# print(ccc)

'''Back Rotor Magnets'''
# revolve top points in circles to the right and left by a certain angle, then extrude that area for the magnet
# Points

bom = motor.copy([(0, 3)]) # Back Outer Magnet Edge Point for Rotation, point tag = 
boml = motor.revolve([bom], 0, 0, 0, 0, 0, 1, sweep_m) # point has tag , curve has tag 
bomr = motor.revolve([boml], 0, 0, 0, 0, 0, -1, 2 * sweep_m) # curve has tag 

bim = motor.copy([(0, 7)]) # Back Inner Magnet Edge Point for Rotation
biml = motor.revolve([bim], 0, 0, 0, 0, 0, 1, sweep_m) # curve has tag 
bimr = motor.revolve([biml], 0, 0, 0, 0, 0, -1, 2 * sweep_m) # curve has tag 

# Connecting Lines
motor.addCircleArc(bomr[0][1],1,boml[0][1],29)
motor.addLine(boml[0][1],biml[0][1],30)
motor.addCircleArc(biml[0][1],1,bimr[0][1],31)
motor.addLine(bimr[0][1],bomr[0][1],32)

# Curved Loop
motor.addCurveLoop([29, 30, 31, 32],7)

# Plane Surface
motor.addPlaneSurface([7],4)

'''Front Rotor Magnets'''
# Points
fom = motor.copy([(0, 12)]) # Front Outer Magnet Edge Point for Rotation, point tag = 25
foml = motor.revolve([fom], 0, 0, 0, 0, 0, 1, sweep_m) # point has tag 26, curve has tag 25
fomr = motor.revolve([foml], 0, 0, 0, 0, 0, -1, 2 * sweep_m) # point has tag 27, curve has tag 26

fim = motor.copy([(0, 16)]) # Front Outer Magnet Edge Point for Rotation, point tag = 28
fiml = motor.revolve([fim], 0, 0, 0, 0, 0, 1, sweep_m) # point has tag 29, curve has tag 27
fimr = motor.revolve([fiml], 0, 0, 0, 0, 0, -1, 2 * sweep_m) # point has tag 30, curve has tag 28

# Connecting Lines
motor.addCircleArc(fomr[0][1],10,foml[0][1],37)
motor.addLine(foml[0][1],fiml[0][1],38)
motor.addCircleArc(fiml[0][1],10,fimr[0][1],39)
motor.addLine(fimr[0][1],fomr[0][1],40)
# Curved Loop
motor.addCurveLoop([37, 38, 39, 40],8)

# Plane Surface
motor.addPlaneSurface([8],5)

'''Stator Core Center'''
# Points
# motor.copy([(0,24),(0,25),(0,26),(0,27)]) # point indexes are 40, 41, 42, 43

# Connecting Lines (using points 24:27 defining inner edges of stator winding)
motor.addLine(24, 25, 41)
motor.addLine(25, 27, 42)
motor.addLine(27, 26, 43)
motor.addLine(26, 24, 44)

# Curved Loop
motor.addCurveLoop([41, 42, 43, 44], 9)

# Plane Surface
motor.addPlaneSurface([9],6)

'''Stator Core Tooth Tips Back'''
# Points
# motor.copy([(0,40),(0,41),(0,42),(0,43)]) # point indexes are 44, 45, 46, 47 on surface of inner stator core

w = wt / cos(sweep_s/2) # horizontal addition to stator tooth tip geometry

motor.addPoint(Riww + w,R-wt,sm - sl/2 - y2, lc, 40)
motor.addPoint(-Riww - w,R-wt,sm - sl/2 - y2,lc, 41)
motor.addPoint(riww + w,r+wt,sm - sl/2 - y2,lc, 42)
motor.addPoint(-riww - w,r+wt,sm - sl/2 - y2,lc, 43)

motor.addPoint(Riww,R-wt,sm - sl/2,lc,100)
motor.addPoint(-Riww,R-wt,sm - sl/2,lc,101)
motor.addPoint(riww,r+wt,sm - sl/2,lc,102)
motor.addPoint(-riww,r+wt,sm - sl/2,lc,103)

# Connecting Lines
# Replicating surface of inner stator core
# motor.addLine(24, 25, 45)
# motor.addLine(25, 27, 46)
# motor.addLine(27, 26, 47)
# motor.addLine(26, 24, 48)
motor.addLine(100, 101, 45)
motor.addLine(101, 103, 46)
motor.addLine(103, 102, 47)
motor.addLine(102, 100, 48)

# Lines of stator tooth tip surface
motor.addLine(40, 41, 49)
motor.addLine(41, 43, 50)
motor.addLine(43, 42, 51)
motor.addLine(42, 40, 52)

# Connecting stator tip to core
motor.addLine(100,40,53)
motor.addLine(101,41,54)
motor.addLine(102,42,55)
motor.addLine(103,43,56)

# Curved Loops
# Z-plane geometry curve loops
motor.addCurveLoop([45, 46, 47, 48],10)
motor.addCurveLoop([49, 50, 51, 52],11)

# Connecting tooth tip start to end of stator core
motor.addCurveLoop([45, 54, -49, -53],12)
motor.addCurveLoop([46, 56, -50, -54],13)
motor.addCurveLoop([47, 55, -51, -56],14)
motor.addCurveLoop([48, 53, -52, -55],15)

# Plane Surface
motor.addPlaneSurface([10],7)
motor.addPlaneSurface([11],8)
motor.addPlaneSurface([12],9)
motor.addPlaneSurface([13],10)
motor.addPlaneSurface([14],11)
motor.addPlaneSurface([15],12)

# Surface Loop
motor.addSurfaceLoop([6,8,9,10,11,12],1)

'''Stator Core Tooth Tips Front'''
# Points
# Shifting points 24:27 the length of the stator
motor.addPoint(Riww,R-wt,sm + sl/2,lc,44)
motor.addPoint(-Riww,R-wt,sm + sl/2,lc,45)
motor.addPoint(riww,r+wt,sm + sl/2,lc,46)
motor.addPoint(-riww,r+wt,sm + sl/2,lc,47)

# Shifting points 40:43 the length of the stator and depth of stator tooth angled portion
motor.addPoint(Riww + w,R-wt,sm + sl/2 + y2, lc, 48)
motor.addPoint(-Riww - w,R-wt,sm + sl/2 + y2,lc, 49)
motor.addPoint(riww + w,r+wt,sm + sl/2 + y2,lc, 50)
motor.addPoint(-riww - w,r+wt,sm + sl/2 + y2,lc, 51)

# Connecting Lines
# Replicating surface of inner stator core
motor.addLine(44, 45, 57)
motor.addLine(45, 47, 58)
motor.addLine(47, 46, 59)
motor.addLine(46, 44, 60)

# Lines of stator tooth tip surface
motor.addLine(48, 49, 61)
motor.addLine(49, 51, 62)
motor.addLine(51, 50, 63)
motor.addLine(50, 48, 64)

# Connecting stator tip to core
motor.addLine(44,48,65)
motor.addLine(45,49,66)
motor.addLine(46,50,67)
motor.addLine(47,51,68)

# Curved Loops
# Z-plane geometry curve loops
motor.addCurveLoop([57, 58, 59, 60],16)
motor.addCurveLoop([61, 62, 63, 64],17)

# Connecting tooth tip start to end of stator core
motor.addCurveLoop([57, 66, -61, -65],18)
motor.addCurveLoop([58, 68, -62, -66],19)
motor.addCurveLoop([59, 67, -63, -68],20)
motor.addCurveLoop([60, 65, -64, -67],21)

# Plane Surface
motor.addPlaneSurface([16],13)
motor.addPlaneSurface([17],14)
motor.addPlaneSurface([18],15)
motor.addPlaneSurface([19],16)
motor.addPlaneSurface([20],17)
motor.addPlaneSurface([21],18)

# Surface Loop
motor.addSurfaceLoop([13,14,15,16,17,18],2)

'''Extruding/Volumes'''
# Rotor Back Plate
motor.extrude([(2,1)],0,0,-1,[8],[rl]) 

# Rotor Back Plate Magnets
brgeom = [0 for i in range(2*p)] # Back Rotor Geometry copy
bmex = [0 for i in range(2*p)] # Back Rotor Rotation

bmex[0] = motor.extrude([(2,4)],0,0,mt,[8],[1]) # extrusion of magnet at top/center (initial for rotation)

for n in range(1,2*p):
    brgeom[n] = motor.copy([bmex[0][1]])
    bmex[n] = motor.rotate([(brgeom[n])], 0, 0, 0, 0, 0, 1, n * pi / p)

# Rotor Front Plate
motor.extrude([(2,2)],0,0,1,[8],[rl]) 

# Rotor Front Plate Magnets
frgeom = [0 for i in range(2*p)] # Front Rotor Geometry copy
fmex = [0 for i in range(2*p)] # Front Rotor Rotation

fmex[0] = motor.extrude([(2,5)],0,0,-mt,[8],[1]) # extrusion of magnet at top/center (initial for rotation)

for n in range(1,2*p):
    frgeom[n] = motor.copy([fmex[0][1]])
    fmex[n] = motor.rotate([(frgeom[n])], 0, 0, 0, 0, 0, 1, n * pi / p)

# Stator Winding
swgeom = [0 for i in range(t)] # Stator Geometry copy in positive z direction
swex = [0 for i in range(t)] # Stator Rotation in positive z direction

swex[0] = motor.extrude([(2,3)],0,0,sl,[8],[1])

# gmsh.model.setColor([swexp[0][:]],255, 0, 0) # Red
# gmsh.model.setColor([swexn[0][:]],255, 0, 0) # Red

for n in range(1,t):
    swgeom[n] = motor.copy([swex[0][1]])
    swex[n] = motor.rotate([(swgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)

# Stator Core
scgeom = [0 for i in range(t)] # Stator Core Geometry copy
scex = [0 for i in range(t)] # Stator Core Rotation

stbgeom = [0 for i in range(t)] # Stator Tip Back Geometry copy
stbex = [0 for i in range(t)] # Stator Tip Back Rotation
stbvgeom = [0 for i in range(t)] # Stator Tip Back Volume Geometry copy
stbvol = [0 for i in range(t)] # Stator Tip Back Volume Rotation

stfgeom = [0 for i in range(t)] # Stator Tip Front Geometry copy
stfex = [0 for i in range(t)] # Stator Tip Front Rotation
stfvgeom = [0 for i in range(t)] # Stator Tip Front Volume Geometry copy
stfvol = [0 for i in range(t)] # Stator Tip Front Volume Rotation


scex[0] = motor.extrude([(2,6)],0,0,sl,[8],[1])
stbex[0] = motor.extrude([(2,8)],0,0,-y3,[8],[1])
stfex[0] = motor.extrude([(2,14)],0,0,y3,[8],[1])

stbvol[0] = motor.addVolume([1])
stfvol[0] = motor.addVolume([2])


for n in range(1,t):
    scgeom[n] = motor.copy([scex[0][1]])
    scex[n] = motor.rotate([(scgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)

    stbvgeom[n] = motor.copy([(3,stbvol[0])])
    stbvol[n] = motor.rotate([(stbvgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)

    stbgeom[n] = motor.copy([stbex[0][1]])
    stbex[n] = motor.rotate([(stbgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)

    stfvgeom[n] = motor.copy([(3,stfvol[0])])
    stfvol[n] = motor.rotate([(stfvgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)

    stfgeom[n] = motor.copy([stfex[0][1]])
    stfex[n] = motor.rotate([(stfgeom[n])], 0, 0, 0, 0, 0, 1, n * 2 * pi/t)


motor.synchronize()

gmsh.option.setNumber('Mesh.CharacteristicLengthMin', 0.05)    
gmsh.option.setNumber('Mesh.CharacteristicLengthMax', 0.05)

gmsh.model.mesh.generate(3)

gmsh.write("3D_motor_mesh.msh")
# gmsh.write("motortest.msh") # when doing code changes/fixing

# print(stfvol[0])
# print(stfvgeom)

gmsh.finalize()