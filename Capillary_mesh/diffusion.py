import numpy as np
from dolfin import *
from IPython import embed

mesh = Mesh('vasculature_volmax2000_mesh.xml')



vess = MeshFunction('int', mesh, 'vasculature_volmax2000_markers.xml')
vessels = MeshFunction('size_t', mesh, 1)


#vessels.array()[vess.where_equal(1)] = 1

#File('2Dvessels.pvd')<<vessels

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

D = 1e-4

bc = Constant(0)

a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx

A = assemble(a)
b = assemble(L)


u_ = interpolate(Constant(1), V)
bc = Constant(0)

bcs = [DirichletBC(V, bc, vessels, 1)]





