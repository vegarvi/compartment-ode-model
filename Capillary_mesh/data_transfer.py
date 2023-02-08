# from fenics_ii.trace_tools.embedded_mesh import EmbeddedMesh
from dolfin import Mesh, MeshFunction, File, assemble, dx, Constant
from branching import find_branches, Walker, vertex_to_DG0_foo, smooth_data
from collections import Counter
import numpy as np

mesh_coarse = Mesh('vasc_mesh.xml')
X_c = mesh_coarse.coordinates()
terminals_c, branches_c = find_branches(mesh_coarse)

# FIXME: smooth data
data = MeshFunction('int', mesh_coarse, 'widths.xml')
smooth_data(data)

# Conversion to something FEniCS_ii can use
foo = vertex_to_DG0_foo(data)

File('Vessels.pvd')<<foo

S_c = np.pi*assemble(foo*dx())*1e-12   #mum^2 --> mm^2 --> m^2 = 1e-12

V_block = 1.1*0.7*0.7*1e-3  # cm^3
V_brain = 1300              # cm^3

n = V_brain/V_block


S_brain = S_c*n

L = assemble(Constant(1)*dx(domain=mesh_coarse))*1e-6*n


print('Human vasculature surface area: ', S_brain)
print('Human vasculature total length: ', L)


#File('wigths_orig.pvd') << foo

# # Fine 1d space
# mesh3d = Mesh('vasculature_volmax2000_mesh.xml')
# markers = MeshFunction('size_t', mesh3d, 'vasculature_volmax2000_markers.xml')
# mesh_fine = EmbeddedMesh(mesh3d, markers, 1).mesh
# X_f = mesh_fine.coordinates()
# terminals_f, branches_f = find_branches(mesh_fine)

# # I want to prevent
# assert all(tuple(t[::-1]) not in terminals_c for t in terminals_c)
# assert all(tuple(t[::-1]) not in terminals_f for t in terminals_f)

# # The first assumption is that we have same number of branches
# assert len(branches_c) == len(branches_f)

# # What about bounding boxes
# assert np.linalg.norm(X_c.min(0) - X_f.min(0)) < 1E-8
# assert np.linalg.norm(X_c.max(0) - X_f.max(0)) < 1E-8

# # There should be a collision of terminals
# all_c_terminals = set(sum(terminals_c, ()))
# all_f_terminals = set(sum(terminals_f, ()))
# assert len(all_c_terminals) == len(all_f_terminals)

# # Map from coarse to fine vertices
# vertex_map = {}
# for c in all_c_terminals:
#     xc = X_c[c]
#     f = min(all_f_terminals, key=lambda f: np.linalg.norm(xc-X_f[f]))
#     assert np.linalg.norm(xc - X_f[f]) < 1E-8, np.linalg.norm(xc - X_f[f]) 
#     vertex_map[c] = f

    
# def collides((c0, c1), (f0, f1)):
#     '''
#     Branch collides as 1 if vertex(c0) == vertex(f0) and likewise for c1, f1
#                     as -1 if vertex(c0) == vertex(f1)
#     '''
#     if vertex_map[c0] == f0 and vertex_map[c1] == f1:
#         return 1
    
#     elif vertex_map[c1] == f0 and vertex_map[c0] == f1:
#         return -1

#     else:
#         return 0

# # Establish correspondence between branches. This is quite funny problem
# # as the 2 end points DO NOT necessary characteize the branch. So some
# # extra work is required.
# ct_count = Counter(terminals_c)
# ft_count = Counter(terminals_f)

# wc = Walker(mesh_coarse)
# wf = Walker(mesh_fine)

# ordered_branches_f = []
# unpaired_fbranches = set(range(len(terminals_f)))
# for cindex, tc in enumerate(terminals_c):
#     # Easy
#     if ct_count[tc] == 1:
#         found = 0
#         for index in unpaired_fbranches:
#             tf, bf = terminals_f[index], branches_f[index]
#             found = collides(tc, tf)
        
#             if found == 1:
#                 ordered_branches_f.append(bf)
#                 unpaired_fbranches.remove(index)
#                 break
        
#             if found == -1:
#                 ordered_branches_f.append(list(reversed(bf)))
#                 unpaired_fbranches.remove(index)
#                 break
#         assert found
#         # Now tf fount
#         assert ft_count[tf] == 1
#     # Loop like case
#     else:
#         bc = branches_c[cindex]
#         # Now we want to look for match
#         for _, d_c in wc.walk(bc): continue

#         # Look for all the coliding branches; wishful thingkin is that
#         # there is only one branch which matches no distance
#         for index in unpaired_fbranches:
#             tf, bf = terminals_f[index], branches_f[index]
#             found = collides(tc, tf)
        
#             if found == 0: continue

#             # Step hoping to match on distance
#             for _, d_f in wf.walk(bf): continue

#             if abs(d_c - d_f) < 1E-4:
#                 ordered_branches_f.append(bf if found == 1 else list(reversed(bf)))
#                 unpaired_fbranches.remove(index)
#                 break
# # we found everybody?
# assert not unpaired_fbranches
                    
# branches_f = ordered_branches_f

# # Okay at this point we have reduced the problem to extending the data
# # to branches; much smaller smaller problem (divide). Now conquer

# data_f = VertexFunction('double', mesh_fine, 0)
# for bc, bf in zip(branches_c, branches_f):
#     coarse_nodes = [p[0] for p in wc.walk(bc)]
#     fine_nodes, cum_dist = list(zip(*list(wf.walk(bf))))

#     # Here =, <, >
#     # This the easiest case as we can feed data without no interpolation
#     if len(coarse_nodes) == len(fine_nodes):
#         for c in coarse_nodes:
#             xc = X_c[c]
#             f = min(fine_nodes, key=lambda f: np.linalg.norm(xc-X_f[f]))
#             assert np.linalg.norm(xc - X_f[f]) < 1E-8, np.linalg.norm(xc - X_f[f])
#             # No prolongation needed
#             data_f[int(f)] = float(data[int(c)])
#     else:
#         assert len(coarse_nodes) < len(fine_nodes)
#         fill_nodes = set(fine_nodes)
#         for c in coarse_nodes:
#             xc = X_c[c]
#             f = min(fine_nodes, key=lambda f: np.linalg.norm(xc-X_f[f]))

#             data_f[int(f)] = float(data[int(c)])
#             fill_nodes.remove(f)
            
#         # For the remaining nodes interpolation is needed
#         # Look for for neighbors of each that can DONATE
#         for f in fill_nodes:
#             this = fine_nodes.index(f)
#             # Setup stuff for interpolation
#             left = this
#             for left in reversed(fine_nodes[:this]):
#                 if left not in fill_nodes:
#                     break
#             vLeft = data[int(left)]
#             dLeft = cum_dist[fine_nodes.index(left)]
#             # Setup stuff for interpolation
#             right = this
#             for right in fine_nodes[this+1:]:
#                 if right not in fill_nodes:
#                     break
#             vRight = data[int(right)]
#             dRight = cum_dist[fine_nodes.index(right)]

#             d = cum_dist[this]
#             data_f[int(f)] = vLeft*(dRight-d)/(dRight-dLeft) + vRight*(d-dLeft)/(dRight-dLeft)

# # Finally
# bar = vertex_to_DG0_foo(data_f)
# File('wigths_smooth_fine.pvd') << bar

# # --------------------------------------------------------------------

# import matplotlib.pyplot as plt

# # Plot data on both meshes
# plt.figure()
# # Let's see how badly the data oscillates
# for branch in branches_c:
#     x, y = [], []
#     for p, d in wc.walk(branch):
#         x.append(d)
#         y.append(data[int(p)])
#     plt.plot(x, y)

# plt.figure()
# # After smoothing
# for branch in branches_f:
#     x, y = [], []
#     for p, d in wf.walk(branch):
#         x.append(d)
#         y.append(data_f[int(p)])
#     plt.plot(x, y)
