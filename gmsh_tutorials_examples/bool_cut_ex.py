import gmsh
from math import sin, cos, tan, pi

'''
Showing boolean cut with matching and non-matching edges
'''

theta = pi/10
theta_p = 2 * theta

lc = 1e-2

r = .75 # inner radius
R = 1.2 # outer radius
D = 2.5 # Domain radius

gmsh.initialize() # gmsh must be initialized in Python before using functions
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("Boolean_Cut") # Creating new model called motor
factory = gmsh.model.occ

# Points
factory.addPoint(0,0,0,lc,1)
factory.addPoint(r,0,0,lc,2)
factory.addPoint(r*cos(theta),r*sin(theta),0,lc,3)
factory.addPoint(R,0,0,lc,4)
factory.addPoint(R*cos(theta),R*sin(theta),0,lc,5)
factory.addPoint(D,0,0,lc,6)
# factory.addPoint(D*cos(theta),D*sin(theta),0,lc,7)
factory.addPoint(D*cos(theta_p),D*sin(theta_p),0,lc,7)

# Lines
factory.addLine(2,4,1)
factory.addCircleArc(4,1,5,2)
factory.addLine(5,3,3)
factory.addCircleArc(3,1,2,4)

factory.addLine(1,6,5)
factory.addCircleArc(6,1,7,6)
factory.addLine(7,1,7)

# Curve Loop
factory.addCurveLoop([1,2,3,4],1)
factory.addCurveLoop([5,6,7],2)

# Surfaces
factory.addPlaneSurface([1],1)
factory.addPlaneSurface([2],2)

# Boolean Cut
factory.cut([(2,2)],[(2,1)],3,removeObject=True,removeTool=True)


factory.synchronize()

gmsh.model.mesh.generate(2)

# gmsh.write("boolean_cut.msh")
gmsh.write("boolean_cut.vtk")

gmsh.finalize()