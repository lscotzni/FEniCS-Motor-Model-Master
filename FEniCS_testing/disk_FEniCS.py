import dolfin as do
import numpy as np
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh
from global_var import Ri, Ro_Ri, CCC_Ri, CCR_Ri

do.set_log_level(1)

# SETTING PROBLEM PARAMETERS (Ri is Inner Radius)
Ro      = Ro_Ri * Ri        # Outer Radius
CCC     = CCC_Ri * Ri       # Center of Current Circle
CCR     = CCR_Ri * Ri       # Radius of Current Circle

I       = 100                # Current (A)
W_area  = np.pi*(CCR)**2    # Winding Area (mm^2)

# MESH IMPORT
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh(
    prefix="disk_mesh",
    dim=2,
    subdomains=True
)

# SETTING CLASSES FOR BOUNDARY CONDITIONS
class InnerBoundary(do.SubDomain):
    def inside(self,x,on_boundary):
        # return on_boundary and np.sqrt(x[0]**2 + x[1]**2) < Ri + do.DOLFIN_EPS
        return on_boundary
        # return on_boundary and np.sqrt(x[0]**2 + x[1]**2) > Ro - do.DOLFIN_EPS
        # return boundaries_mf[0], boundaries_mf[1]

class PeriodicBoundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return x[0] > do.DOLFIN_EPS and do.near(x[1],0)

    def map(self,x,y):
        y[0] = -x[0]
        y[1] = x[1]

inner_bound         = InnerBoundary()
periodic_bound      = PeriodicBoundary()

V       = do.FunctionSpace(mesh,'P',1,constrained_domain=periodic_bound)
dx      = do.Measure('dx', domain = mesh, subdomain_data = subdomains_mf)                
A_inner = do.Constant(0)
bc_i    = do.DirichletBC(V,A_inner,inner_bound)

J = do.Constant(I / W_area)

class Permeability(do.UserExpression):
    def __init__(self,subdomains_mf,**kwargs):
        super().__init__(**kwargs)
        self.subdomains_mf = subdomains_mf
    def eval_cell(self,values,x,cell):
        if self.subdomains_mf[cell.index] == 1:
            values[0] = 1.26e-2
            # values[0] = 6.3e-3
        elif self.subdomains_mf[cell.index] == 2:
            values[0] = 6.3e-3
        else:
            values[0] = 1.26e-6

mu      = Permeability(subdomains_mf, degree=1)
# A_z     = do.TrialFunction(V)
A_z     = do.Function(V)
v       = do.TestFunction(V)
a       = (1/mu)*do.dot(do.grad(A_z),do.grad(v))*dx
L       = J*v*dx(3) - J*v*dx(4)
# A_z     = do.Function(V)

do.solve(a-L==0,A_z,bc_i)

asdf = do.assemble(a - L).get_local()

print(do.assemble(a - L).get_local())
print('norm = ',do.norm(do.assemble(a - L)))
print(len(do.assemble(a - L).get_local()))
exit()


W       = do.VectorFunctionSpace(mesh,'P',1)
B       = do.project(do.as_vector((A_z.dx(1),-A_z.dx(0))),W)



''' ================== START WORK FOR SPATIAL PERMEABILITY ================== '''


def spatial_permeability(index,B=0):
    if index == 1:
        permeability_function = 6.3e-3 + B*1e-2
    elif index == 2:
        permeability_function = 6.3e-3
    else: # Constant for copper
        permeability_function = 1.26e-6
    
    return permeability_function

def integrand(permeability_function,solution,test_function):
    # permeability_function represents mu in terms of B ====> mu(B)
    # solution should 
    return (1/permeability_function) * do.dot(do.grad(solution),do.grad(test_function))

print(len(B.vector().get_local()))
# exit()

asdf       = integrand(spatial_permeability(index=1,B=B.vector().get_local()),A_z,v)*dx(1) + integrand(spatial_permeability(index=2),A_z,v)*dx(2) + \
                    integrand(spatial_permeability(index=3),A_z,v)*dx(3) + integrand(spatial_permeability(index=3),A_z,v)*dx(4)

''' ================== END WORK FOR SPATIAL PERMEABILITY ================== '''

do.plot(B,linewidth=20)

vtkfile_A_z = do.File('solutions/Magnetic_Vector_Potential.pvd')
vtkfile_B = do.File('solutions/Magnetic_Flux_Density.pvd')
vtkfile_A_z << A_z
vtkfile_B << B

plt.show()