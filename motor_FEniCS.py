from fenics import *
import dolfin as do
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh_from_xdmf

# Import parameters from mesh
# from motor_mesh_2D import p,s

set_log_level(1)

# Stator winding domains have been defined; these can have specific J values
# relating to current density in Maxwell's equations; 3-phase means phase shift of pi/3 rad
# between each current
# Boundary conditions set on x-axis and edge of domain
# individual domains have permeability depending on area (air, etc)

# We have listed subdomains in the GMSH file depending on position, indexed from 0 to 10
# Air gap: 0; Rotor core: 1; Stator core: 2; Odd Magnets: 3; Even Magnets: 4
# Stator Windings: A: 5,6; B: 7,8; C: 9,10

tol = 1e-15
print(DOLFIN_EPS)

# Setting geometrical conditions for boundaries
class Domain_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[1] > DOLFIN_EPS

# Setting geometry for periodic boundary condition
class Semicircle_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] > -DOLFIN_EPS and near(x[1],0)
        # Taking right side of semicircle boundary as original "target domain"

    def map(self,x,y): # Mapping to new boundary, w/ same conditions for y-axis and negative 
        y[0] = -x[0]
        y[1] = x[1]

domain_boundary = Domain_Boundary() # Far-Field Boundary Condition (approaches zero)
pbc = Semicircle_Boundary() # Periodic Boundary Conditions
# -------------------------------------------------------------------------------------------
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh_from_xdmf(
    prefix="motor_mesh",
    dim=2,
    subdomains=True
)


# -------------------------------------------------------------------------------------------
# mesh = do.Mesh()
# # with do.XDMFFile("2d_motor.xdmf") as infile:
# with do.XDMFFile("2d_motor_terminal.xdmf") as infile:
#     infile.read(mesh)
# mvc = MeshValueCollection("size_t",mesh,2)
# with do.XDMFFile("2d_motor_mf.xdmf") as infile:
#     infile.read(mvc, "name_to_read")
# mf = cpp.mesh.MeshFunctionSizet(mesh,mvc)

# domains = MeshFunction("size_t", mesh, mesh.topology().dim())
# dx = Measure("dx",domain=mesh,subdomain_data=mf)
# print(assemble(Constant(1)*dx))
# print(assemble(Constant(1)*dx(1)))

# -------------------------------------------------------------------------------------------
# mesh = do.Mesh("2d_motor.xml")
# mvc = MeshValueCollection("size_t",mesh,2)
# mf = cpp.mesh.MeshFunctionSizet(mesh,mvc)
    
# subdomains = do.MeshFunction("size_t", mesh, mesh.topology().dim())
# # dx = Measure("dx",domain=mesh, subdomain_data=mf)   
# dx = Measure('dx',subdomain_data=mf)
# # print(assemble(Constant(1)*dx(domain=mesh)))
# # print(assemble(Constant(1)*dx(1)))
# dx0 = Measure('dx',domain = mesh, subdomain_data = subdomains,subdomain_id = 0) # Airgap integration measure
# dx1 = Measure('dx',domain = mesh, subdomain_data = subdomains,subdomain_id = 2) # Airgap integration measure
# print(assemble(Constant(1)*dx0))
# print(assemble(Constant(1)*dx1))

# boundaries = do.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# boundaries.set_all(0)
# domain_boundary.mark(boundaries,1)



# fig = plt.figure(1,figsize = (9,8))
# ax = fig.gca(projection="3d")
# ax.view_init(elev=90., azim=-90.)
# plt.axis([-4.5,4.5,0,4.25])
do.plot(mesh)
plt.show()