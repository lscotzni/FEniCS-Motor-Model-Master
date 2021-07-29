from fenics import *
import dolfin as do
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh_from_xdmf

# Import point data from mesh
# from rectangle_mesh import l, w   # l is the height (y) and w is width (x)
l = 2
w = 6

'''
This FEniCS Example shows a simple rectangular mesh imported 
from gmsh with Dirichlet BC along the x-axis and insulated conditions 
on the y-axis boundaries

Steady, 2D, homogeneous heat equation: Laplace(u) = 0
'''
set_log_level(1)

# Defining Boundary/Domain Classes
tol = 1e-15

class Top_Temp(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],2*l)

class Bottom_Temp(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],0)

class Left_Temp(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] < tol

class Right_Temp(do.SubDomain):
    def inside(self,x,on_boundary): 
        return on_boundary and x[0] > w - tol

class Mid_Temp(do.SubDomain):
    def inside(self,x,on_boundary):
        return near(x[1],l)

class Upper_Surface(do.SubDomain):
    def inside(self,x,on_boundary):
        return x[1] > l - tol

class Lower_Surface(do.SubDomain):
    def inside(self,x,on_boundary):
        return x[1] < l + tol

left = Left_Temp()
right = Right_Temp()
top = Top_Temp()
bottom = Bottom_Temp()
mid = Mid_Temp()
upper = Upper_Surface()
lower = Lower_Surface()

# mesh,boundaries_mf, subdomains_mf, association_table = import_mesh_from_xdmf(
#     prefix = "rect_mesh",
#     dim = 2,
#     subdomains = True
# )

# # Assigning Mesh Domains/Boundaries
# mesh = do.Mesh('rect_mesh.xml')
# # mesh = RectangleMesh(Point(0.0, 0.0), Point(w, 2*l), 10, 10)

# # -------------------------------------------------------------------------------------------

domains = do.MeshFunction("size_t", mesh, mesh.topology().dim())
boundaries = do.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
# third argument gives the topological dimension of the MeshFunction, 
# which is the topological dimension of our mesh minus 1.
# mesh.topology().dim() = 2; size_t represents unsigned integer (uint)
'''Domains'''
domains.set_all(0) # Marking all domain in mesh as 0
upper.mark(domains,1)
lower.mark(domains,2)

'''Boundaries'''
boundaries.set_all(0)
left.mark(boundaries,1)
right.mark(boundaries,2)
top.mark(boundaries,3)
bottom.mark(boundaries,4)
mid.mark(boundaries,5)

dx = Measure('dx', domain=mesh, subdomain_data=boundaries)

# class Material(UserExpression):
#     def __init__(self,domains,m0,m1,**kwargs):
#         self.domains = domains
#         self.m0 = m0
#         self.m1 = m1

# def eval(self,value,x):
#     if self.materials[cell.index] == 1:
#         value[0] = self.m0
#     elif self.materials[cell.index] == 2:
#         value[0] = self.m1

# m0, m1 = 5, 25

# m = Material(domains, m0, m1, degree=0)

# Defining function space and basis functions
V = FunctionSpace(mesh,'P',1) # 'P' is the same as CG & Lagrange
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Boundary Conditions
L_T = Constant(0.0) # Prescribed Left Side Temperature (can change to an expression or function)
R_T = Constant(3.0) # Prescribed Right Side Temperature
# T_T = Expression("x[0]/2", degree = 1)
# B_T = Expression("x[0]/2", degree = 1)
T_T = Expression("x[0]*(x[0]-5)/2", degree = 2)
B_T = Expression("-x[0]*(x[0]-7)/2", degree = 2)
M_T = Expression("x[0]*x[0]/12", degree = 2)

bcs = [
    DirichletBC(V,L_T,left),
    DirichletBC(V,R_T,right),
    DirichletBC(V,T_T,top),
    DirichletBC(V,B_T,bottom),
    DirichletBC(V,M_T,mid)
]

# Compute Solution
u = Function(V)
solve(a == L, u, bcs)

fig = plt.figure(1,figsize = (9,8))
ax = fig.gca(projection="3d")
ax.set_aspect('auto')
ax.view_init(elev=90., azim=-90.)
plt.axis([-1,w+1,-1,2*l+1])

plot(u)
plot(mesh)
plt.show()