from dolfin import *
import numpy as plt

from IPython import embed

mesh = Mesh('2Didealized.xml')          # MESH IN mm
mesh = refine(mesh)
mesh = refine(mesh)

L = 1.12
H = 1.12
d = 0.28
R = 0.015


points = []
ART = [0,2,5,7,8,10,13,15]
VEI = [1,3,4,6,9,11,12,14]
for i in range(1,5):
    for j in range(1,5):
        points.append([i*d - d/2, j*d-d/2])

artr = [points[i] for i in ART]
vein = [points[i] for i in VEI]


class walls(SubDomain):
    def inside(self,x,on_bnd):
        return (x[0] < DOLFIN_EPS or x[0] > L - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > H - DOLFIN_EPS) and on_bnd

class arteries(SubDomain):
    def inside(self,x,on_bnd):
        for p in artr:
            r = plt.sqrt((x[0] - p[0])**2 + (x[1] - p[1])**2)
            if (r < 2*R) and on_bnd:
                return True
        return False

class veins(SubDomain):
    def inside(self,x,on_bnd):
        for p in vein:
            r = plt.sqrt((x[0] - p[0])**2 + (x[1] - p[1])**2)
            if (r < 2*R) and on_bnd:
                return True
        return False


bnd = MeshFunction('size_t', mesh, 1)


walls().mark(bnd,1)
arteries().mark(bnd,2)
veins().mark(bnd,3)

File('bnd.pvd')<<bnd

V = FunctionSpace(mesh, 'CG', 1)

p = TrialFunction(V)
q = TestFunction(V)

K = 2e-17*1e6       # Permeability      [mm^2]
mu = 0.7*1e-3       # Viscosity         [Pa*s] 
    
p_a = Constant(50*133.33)       # Arteriole pressure [Pa]
p_v = Constant(15*133.33)       # Venule pressure    [Pa]

ds = Measure('ds')(subdomain_data=bnd)
bca = DirichletBC(V,p_a,bnd,2)   	    # Arteriole bc
bcv = DirichletBC(V,p_v,bnd,3)          # Venule bc
bcs = [bca, bcv]



A = assemble(K/mu*inner(grad(p),grad(q))*dx)
b = assemble(Constant(0)*q*dx)
[bc.apply(A,b) for bc in bcs]
p_ = Function(V)

solve(A, p_.vector(), b)
n = FacetNormal(mesh)
flow_cube = assemble(K/mu*dot(grad(p_),n)*ds(2))     # [mm^2/s] ( mm^3/s per mm depth)

area = L*H                                 # Rectangle area                       [mm^2]
brain_vol = 1300*1e3                       # Brain volume                         [mm^3]
ratio = brain_vol/area                     # Length/number of cubes in full brain [mm]

Q = ratio*flow_cube*1e-9                   # total flow in m^3/s (from mm^3/s)

print(Q)

R = Q**-1*(float(p_a) - float(p_v))        # Resistance in Pa*s/m^3
p_convert = 1e6*133*60.                    # converting from mmHg/(ml/min) to Pa *s /m^3

print(R/p_convert)                         # Resistance in mmHg/(mL/min)


File('results_darcy/pressure.pvd')<< p_













