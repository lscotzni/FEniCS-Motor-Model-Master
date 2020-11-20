import gmsh
import numpy as np 

gmsh.initialize()
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("asdf") # Creating new model called motor
asdf = gmsh.model.occ

lc = 1e-2

a = 2

asdf.addPoint(0,0,0,lc,1)
asdf.addPoint(a,0,0,lc,2)
asdf.addPoint(a,a,0,lc,3)
asdf.addPoint(0,a,0,lc,4)
asdf.addPoint(0,a,0,lc,5)

# gmsh.model.mesh.removeDuplicateNodes()

asdf.addLine(1,2,1)
asdf.addLine(2,3,2)
asdf.addLine(3,4,3)
asdf.addLine(4,1,4)
asdf.addLine(5,1,5)

gmsh.model.occ.removeAllDuplicates()

asdf.addCurveLoop([1, 2, 3, 4],1)
asdf.addCurveLoop([-4, -3, -2, -1],2)

s1 = asdf.addPlaneSurface([1],1)
s2 = asdf.addPlaneSurface([2],2)

gmsh.model.occ.removeAllDuplicates()

asdf.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("duplicate_test.msh")

gmsh.finalize()