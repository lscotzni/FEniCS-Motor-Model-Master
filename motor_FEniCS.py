from fenics import *
import dolfin as do
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh
from math import sin, cos, pi

set_log_level(1)

# Stator winding domains have been defined; these can have specific J values
# relating to current density in Maxwell's equations; 3-phase means phase shift of 2pi/3 rad
# between each current
# Boundary conditions set on x-axis and edge of domain
# individual domains have permeability depending on area (air, etc)

tol = 1e-15

s = 9
p = 8


theta = 0
# Imax = 282.8 # (A) Max current in wires defined in paper
Imax = 0 # (A) No Load Current
# Nstp = 39 # Number of turns per phase
Nstp = 1 

theta_t = pi/s
theta_sso = .5 * theta_t
Rsy = 103.
Rssi = 83.

theta_p = 2*pi/2/p # Angular sweep of each Rotor Slice
theta_m = .78 * theta_p # Angular sweep of magnet
theta_b = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
theta_g = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet side piece

Rtm = 79.
Rbm = 74.5

Hc = 800 # Coercivity of Neodymium magnet in kA/m (ranges between 800 and 950)
lm = Rtm - Rbm # Thickness of magnet (mm)

winding_area = theta_sso/2/2 * (Rsy**2 - Rssi**2) * 10**(-6)
magnet_area = (Rtm**2 - Rbm**2) * theta_m/2

# Setting geometrical conditions for boundaries
class Domain_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[1] > DOLFIN_EPS

# Setting geometry for periodic boundary condition
class Periodic_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and near(x[1],0.)
        # Taking right side of semicircle boundary as original "target domain"

    def map(self,x,y): # Mapping to new boundary, w/ same conditions for negative part of x-axis  
        y[0] = -x[0]
        y[1] = x[1]

domain_bound = Domain_Boundary() # Far-Field Boundary Condition (approaches zero)
semicircle_bound = Periodic_Boundary() # Periodic Boundary Conditions
# -------------------------------------------------------------------------------------------
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh(
    prefix="motor_mesh",
    dim=2,
    subdomains=True
)
# boundaries_mf: Mesh Function object for boundaries
# subdomains_mf: Mesh Function object for subdomains

V = FunctionSpace(mesh,'P',1, constrained_domain = Periodic_Boundary())

dx = Measure('dx', domain = mesh, subdomain_data = subdomains_mf)
dS = Measure('dS', domain = mesh, subdomain_data = boundaries_mf)


'''
----- Domain list for dx: -----
Rotor Core: 1
Stator Core: 2
Outer Motor Domain: 3
Inner Motor Domain: 4
Middle Air Gap: 5
Right Air Slot Domains (next to magnet): 6 to 13
Left Air Slot Domains (next to magnet): 14 to 21
Magnet Domains: 22 to 29
Right Stator Winding Domains: 30 to 38
Left Stator Winding Domains: 39 to 47

----- Boundary list for ds: -----
Magnet Edges (for line current density): 100 to 115
'''

# Boundary Conditions
A_FF = Constant(0) # Far Field Magnetic Vector Potential (decays to zero)
# Need to figure out how to assign the periodic BC above to the x-axis 
# (may be done w/ constrained domain in Function space definition)
bc0 = DirichletBC(V,A_FF,Domain_Boundary())

# Current Densities and Winding Indices
J1_value = Imax * Nstp * cos(theta) / winding_area # Current Density 1 value
J2_value = Imax * Nstp * cos(theta + 2*pi/3) / winding_area # Current Density 2 value
J3_value = Imax * Nstp * cos(theta - 2*pi/3) / winding_area # Current Density 3 value
IM_value = Hc * lm # Effective current in magnet (A)
# JM_value = IM_value / magnet_area # Effective current density in magnet (A/m^2)
JM_value = IM_value

J1 = Constant(J1_value) # Current Density 1
J2 = Constant(J2_value) # Current Density 2
J3 = Constant(J3_value) # Current Density 3
JM = Constant(JM_value) # Equivalent Line Current Density of Magnets

WA_N = [30, 40, 32] # Winding Domain indices for current A North
WA_S = [39, 31, 41] # Winding Domain indices for current A South
WB_N = [33, 43, 35] # Winding Domain indices for current B North
WB_S = [42, 34, 44] # Winding Domain indices for current B South
WC_N = [36, 46, 38] # Winding Domain indices for current C North
WC_S = [45, 37, 47] # Winding Domain indices for current C South

MI = [0 for i in range(p)] # Magnet Domain segment indices for inward current (into page, denoted as negative "-")
MO = [0 for i in range(p)] # Magnet Domain segment indices for outward current (out of page, denoted as positive "+")
init_segment = association_table["current density inward line segment 1"] # First segment of magnet from dictionary

for i in range(p):
    if i%2 == 0:
        MI[i] = init_segment + 2*i
        MO[i] = init_segment + 2*i + 1
    else:
        MO[i] = init_segment + 2*i
        MI[i] = init_segment + 2*i + 1

# Define Magnetic Permeability
class Permeability(UserExpression):
    def __init__(self,subdomains_mf,**kwargs):
        super().__init__(**kwargs)
        self.subdomains_mf = subdomains_mf
    def eval_cell(self, values, x, cell):
        if self.subdomains_mf[cell.index] == 1:
            # values[0] == 1 # Insert permeability value for rotor core material
            # values[0] = 1e-5
            values[0] = 6.3e-3
        elif self.subdomains_mf[cell.index] == 2:
            # values[0] == 1 # Insert permeability value for stator core material
            # values[0] = 1e-5
            values[0] = 6.3e-3
        elif (self.subdomains_mf[cell.index] >= 3) and (self.subdomains_mf[cell.index] <= 21): # subdomains_mf values for the air gaps (may remove 4 bc it is the shaft) (ALSO NEED TO INCLUDE AIR SLOTS AT MAGNET)
            # values[0] == 1 # Insert permeability value for air
            # values[0] = 4*pi*1e-7
            values[0] = 6.3e-3
        elif (self.subdomains_mf[cell.index] >= 22) and (self.subdomains_mf[cell.index] <= 29):
            # values[0] = 1.05 * 4*pi*1e-7 # Insert permeability value for magnets (this is defined as mu_rm = 1.05)
            values[0] = 6.3e-3
        elif self.subdomains_mf[cell.index] >= 30:
            # values[0] == 1 # Insert permeability value for wires (copper)
            # values[0] = 1.26e-6
            values[0] = 6.3e-3

mu = Permeability(subdomains_mf, degree=1)

# Define Variational Problem
A_z = TrialFunction(V)
v = TestFunction(V)
a = (1 / mu)*dot(grad(A_z), grad(v))*dx

# Current Densities for each winding phase
JA_N = sum(J1*v*dx(WA_N[i]) for i in range(3))
JA_S = -sum(J1*v*dx(WA_S[i]) for i in range(3))
JB_N = sum(J2*v*dx(WB_N[i]) for i in range(3))
JB_S = -sum(J2*v*dx(WB_S[i]) for i in range(3))
JC_N = sum(J3*v*dx(WC_N[i]) for i in range(3))
JC_S = -sum(J3*v*dx(WC_S[i]) for i in range(3))

# Current Densities for magnets
JM_I = -sum(JM*v('+')*dS(MI[i]) for i in range(p))
JM_O = sum(JM*v('+')*dS(MO[i]) for i in range(p))

L = JA_N + JA_S + JB_N + JB_S + JC_N + JC_S + JM_O + JM_I
# L = JA_N + JA_S + JB_N + JB_S + JC_N + JC_S
# L = Constant(0) * v * dx

''' ------------------------ ANALYSIS METHODS ------------------------ '''
# Solve Variational Problem
A_z = Function(V)


method = 'point_source' # ANALYSIS METHOD

if method == 'base': # ------------------------ Base Method ------------------------
    solve(a == L, A_z, bc0)

elif method == 'point_source': # ------------------------ Point Source Method ------------------------
    LHS, RHS = assemble_system(a,L,bc0)

    theta_lower = theta_b/2 + theta_g # Angular position of lower angular displacement of magnet
    # theta_lower = (theta_b/2 + theta_g) * .5
    Rma = (Rbm + Rtm)/2 # Average (midpoint) of top and bottom magnet surfaces

    for i in range(p):
        print("Iteration: ",i+1)
        theta_higher = theta_lower + theta_m

        print(theta_lower * 180/pi)
        print(theta_higher * 180/pi)
        if i%2 == 0:
            x_in, y_in = Rma * cos(theta_lower), Rma * sin(theta_lower)
            x_out, y_out = Rma * cos(theta_higher), Rma * sin(theta_higher)
        else:
            x_in, y_in = Rma * cos(theta_higher), Rma * sin(theta_higher)
            x_out, y_out = Rma * cos(theta_lower), Rma * sin(theta_lower)

        print(x_in,y_in)
        print(x_out,y_out)

        delta_in = PointSource(V,Point(x_in,y_in),-JM_value)
        delta_in.apply(RHS)

        delta_out = PointSource(V,Point(x_out,y_out),JM_value)
        delta_out.apply(RHS)

        theta_lower += pi/8

    random_PS1 = PointSource(V,Point(0,55),JM_value)
    random_PS1.apply(RHS)

    random_PS2 = PointSource(V,Point(0,65),JM_value)
    random_PS2.apply(RHS)

    solve(LHS,A_z.vector(),RHS)

elif method == 'int_facet': # ------------------------ Interior Facet Method ------------------------
    pass


# Computing magnetic field (B = curl(A))
W = VectorFunctionSpace(mesh,'P',1)
B = project(as_vector((A_z.dx(1), -A_z.dx(0))),W)

'''MESH PLOT'''
# fig = plt.figure(1,figsize = (9,8))
# ax = fig.gca(projection="3d")
# ax.view_init(elev=90., azim=-90.)
# plt.axis([-4.5,4.5,0,4.25])
# do.plot(subdomains_mf)
# do.plot(A_z)
# do.plot(mesh, linewidth = .1)
do.plot(B, linewidth = 20)


vtkfile_A_z = File('solutions/Magnetic_Vector_Potential.pvd')
vtkfile_B = File('solutions/Magnetic_Flux_Density.pvd')
vtkfile_A_z << A_z
vtkfile_B << B

plt.show()

# -------------------------------------------------------------------------------------------
'''TO-DO'''
# Set proper permeability values (and figure out the thing w vacuum)
# Assign the whole boundary at the x-axis 

# Figure out what the debug code means in output:
# DEBUG: [at /Users/lucascotzniovsky/opt/anaconda3/include/dolfin/mesh/MeshFunction.h:485 in operator=()]
# DEBUG: Mesh value collection does not contain all values for all entities

# for dx, we can call dx(n), where n is the index of the domain number in the Mesh Function
# we need this for the stator windings to assign the values of J1,J2,J3 to each winding
# This has to be done manually, but we can arrange a loop for currents in and out of the 