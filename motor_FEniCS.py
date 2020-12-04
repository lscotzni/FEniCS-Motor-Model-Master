from fenics import *
import dolfin as do
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh_from_xdmf
from math import cos, pi

# Import parameters from mesh
# from motor_mesh_2D import p,s

# set_log_level(1)

# Stator winding domains have been defined; these can have specific J values
# relating to current density in Maxwell's equations; 3-phase means phase shift of pi/3 rad
# between each current
# Boundary conditions set on x-axis and edge of domain
# individual domains have permeability depending on area (air, etc)

tol = 1e-15
print(DOLFIN_EPS)

theta = 0
Imax = 282.8 # (A) Max current in wires defined in paper
Nstp = 39 # Number of turns per phase

# Setting geometrical conditions for boundaries
class Domain_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[1] > DOLFIN_EPS

# Setting geometry for periodic boundary condition
class Periodic_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and near(x[1],0)
        # Taking right side of semicircle boundary as original "target domain"

    def map(self,x,y): # Mapping to new boundary, w/ same conditions for negative part of x-axis  
        y[0] = -x[0]
        y[1] = x[1]

domain_bound = Domain_Boundary() # Far-Field Boundary Condition (approaches zero)
semicircle_bound = Periodic_Boundary() # Periodic Boundary Conditions
# -------------------------------------------------------------------------------------------
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh_from_xdmf(
    prefix="motor_mesh",
    dim=2,
    subdomains=True
)
# boundaries_mf: Mesh Function object for boundaries
# subdomains_mf: Mesh Function object for subdomains

V = FunctionSpace(mesh,'P',1, constrained_domain = Periodic_Boundary())


dx = Measure('dx', domain = mesh, subdomain_data = subdomains_mf)
'''
Domain list for dx:
Rotor Core: 1
Stator Core: 2
Outer Motor Domain: 3
Inner Motor Domain: 4
Middle Air Gap: 5
Right Air Slot Domains (next to magnet): 6 to 13
Left Air Slot Domains (next to magnet): 14 to 21
Magnet Domains (next to magnet): 22 to 29
Right Stator Winding Domains: 30 to 38
Left Stator Winding Domains: 39 to 47
'''
# print(assemble(Constant(1)*dx))
# ds = Measure('ds', domain = mesh, subdomain_data = boundaries_mf) # need to confirm, but this is setting measure for the boundaries


# Boundary Conditions
A_FF = Constant(0) # Far Field Magnetic Vector Potential (decays to zero)
# Need to figure out how to assign the periodic BC above to the x-axis 
# (may be done w/ constrained domain in Function space definition)
bc0 = DirichletBC(V,A_FF,Domain_Boundary())

# Current Densities and Winding Indices
J1 = Constant(Imax * Nstp * cos(theta)) # Current Density 1
J2 = Constant(Imax * Nstp * cos(theta + pi/3)) # Current Density 2
J3 = Constant(Imax * Nstp * cos(theta + 2*pi/3)) # Current Density 3

WA_N = [30, 40, 32] # Winding Domain indices for current A North
WA_S = [39, 31, 41] # Winding Domain indices for current A South
WB_N = [33, 43, 35] # Winding Domain indices for current B North
WB_S = [42, 34, 44] # Winding Domain indices for current B South
WC_N = [36, 46, 38] # Winding Domain indices for current C North
WC_S = [45, 37, 47] # Winding Domain indices for current C South

# Define Magnetic Permeability
class Permeability(UserExpression):
    def __init__(self,subdomains_mf,**kwargs):
        super().__init__(**kwargs)
        self.subdomains_mf = subdomains_mf
    def eval_cell(self, values, x, cell):
        if self.subdomains_mf[cell.index] == 1:
            # values[0] == 1 # Insert permeability value for rotor core material
            values[0] = 1e-5
        elif self.subdomains_mf[cell.index] == 2:
            # values[0] == 1 # Insert permeability value for stator core material
            values[0] = 1e-5
        elif self.subdomains_mf[cell.index] >= 3 and self.subdomains_mf[cell.index] <= 21: # subdomains_mf values for the air gaps (may remove 4 bc it is the shaft) (ALSO NEED TO INCLUDE AIR SLOTS AT MAGNET)
            # values[0] == 1 # Insert permeability value for air
            values[0] == 4*pi*1e-7
        elif self.subdomains_mf[cell.index] >= 22 and self.subdomains_mf[cell.index] <= 29:
            values[0] == 1 # Insert permeability value for magnets (this is defined as mu_rm = 1.05)
        else:
            # values[0] == 1 # Insert permeability value for wires (copper)
            values[0] == 1.26e-6

mu = Permeability(subdomains_mf, degree=1)

# Define Variational Problem
A_z = TrialFunction(V)
v = TestFunction(V)
a = (1 / mu)*dot(grad(A_z), grad(v))*dx

JA_N = sum(J1*v*dx(WA_N[i]) for i in range(3))
JA_S = -sum(J1*v*dx(WA_S[i]) for i in range(3))
JB_N = sum(J2*v*dx(WB_N[i]) for i in range(3))
JB_S = -sum(J2*v*dx(WB_S[i]) for i in range(3))
JC_N = sum(J3*v*dx(WC_N[i]) for i in range(3))
JC_S = -sum(J3*v*dx(WC_S[i]) for i in range(3))

L = JA_N + JA_S + JB_N + JB_S + JC_N + JC_S # (sum up these 6 lines above)

# Solve Variational Problem
A_z = Function(V)
solve(a == L, A_z, bc0)

# Computing magnetic field (B = curl(A))
W = VectorFunctionSpace(mesh,'P',1)
B = project(as_vector((A_z.dx(1), -A_z.dx(0))),W)

'''MESH PLOT'''
# fig = plt.figure(1,figsize = (9,8))
# ax = fig.gca(projection="3d")
# ax.view_init(elev=90., azim=-90.)
# plt.axis([-4.5,4.5,0,4.25])
# do.plot(subdomains_mf)
do.plot(A_z)
# do.plot(mesh)
plt.show()

# -------------------------------------------------------------------------------------------
'''TO-DO'''
# Set proper permeability values (and figure out the thing w vacuum)
# Assign the whole boundary at the x-axis 

# Figure out what the debug code means in output:
# DEBUG: [at /Users/lucascotzniovsky/opt/anaconda3/include/dolfin/mesh/MeshFunction.h:485 in operator=()]
# DEBUG: Mesh value collection does not contain all values for all entities

# Error in plotting A_z -> ValueError: z array must not contain non-finite values within the triangulation
# Plot of B shows nothing


# for dx, we can call dx(n), where n is the index of the domain number in the Mesh Function
# we need this for the stator windings to assign the values of J1,J2,J3 to each winding
# This has to be done manually, but we can arrange a loop for currents in and out of the 