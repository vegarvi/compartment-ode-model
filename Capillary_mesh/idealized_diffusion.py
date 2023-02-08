from writereal import exclude_boarder_vessels
from dolfin import *
import numpy as plt

from IPython import embed

mesh = Mesh('real_output_geo_large_veins.xml')          # MESH IN mm
#mesh = refine(mesh)
#mesh = refine(mesh)


x0,y0,l,h = plt.loadtxt('origin_L_H.txt')
R_a = 0.015        # [mm]
R_v = 0.020        # [mm]


artr = exclude_boarder_vessels(plt.loadtxt('arteries.txt'), R_a)
vein = exclude_boarder_vessels(plt.loadtxt('veins.txt'), R_a)


class walls(SubDomain):
    def inside(self,x,on_bnd):
        return (x[0] < DOLFIN_EPS or x[0] > l - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > h - DOLFIN_EPS) and on_bnd

class arteries(SubDomain):
    def inside(self,x,on_bnd):
        for p in artr:
            r = plt.sqrt((x[0] - p[0])**2 + (x[1] - p[1])**2)
            if (r < 2*R_a) and on_bnd:
                return True
        return False

class veins(SubDomain):
    def inside(self,x,on_bnd):
        for p in vein:
            r = plt.sqrt((x[0] - p[0])**2 + (x[1] - p[1])**2)
            if (r < 2*R_v) and on_bnd:
                return True
        return False


bnd = MeshFunction('size_t', mesh, 1)


walls().mark(bnd,1)
arteries().mark(bnd,2)
veins().mark(bnd,3)

File('bnd.pvd')<<bnd


V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

D = Constant(0.7e-11*1e6)      # Diffusion coefficient      [mm^2]  (can be 1e-4 - 1e-3 nm^2)

    
c_v = Constant(0)
c_v2 = Constant(1)

ds = Measure('ds')(subdomain_data=bnd)
bca = DirichletBC(V,c_v,bnd,2)   	    # Arteriole bc
bcv = DirichletBC(V,c_v2,bnd,3)          # Venule bc
bcs = [bca, bcv]


u1 = interpolate(Constant(1), V)

dt = Constant(100)
T = 3600*5


a = inner(u, v)*dx + dt*D*inner(grad(u), grad(v))*dx
L = inner(u1, v)*dx



A = assemble(a)



u_ = Function(V)
t = float(dt)
outfile = File('diffusion/concentration.pvd')

dn = 0
#V_brain = 1300*1e3
#V_csf = 130*1e3
#Area = l*h
#nb0 = assemble(u1*dx)*V_brain/Area
#V_brain = 1300*1e3
#V_csf = 130*1e3
#c_csf0 = 0.8
#Q = 55.0                               # Outflow rate

C0 = assemble(u1*dx)
C1 = C0

while C1/C0 > 0.5:
    #c_csf = -dn/V_csf + c_csf0 - Q*c_csf0*float(dt)/V_csf
    #if t > 300:
    #    dt = Constant(50)
    #print('csf conc: ', c_csf)

    b = assemble(L)
    [bc.apply(A,b) for bc in bcs]
    solve(A, u_.vector(), b)
    #nb1 = assemble(u_*dx)*V_brain/Area
    C1 = assemble(u_*dx)
    
    #dn = float(nb1) - float(nb0)
    u1.assign(u_)
    outfile << u_
    t += float(dt)
    #nb0 = nb1
    #c_csf0 = c_csf
    print('t = ', t, 'Relative concentration: ', C1/C0)

print('Half-life: ', t)
exit()
area = L*H                                           # Rectangle area                       [mm^2]
brain_vol = 1300*1e3                                 # Brain volume                         [mm^3]
ratio = brain_vol/area*1e-3                          # Length/number of cubes in full brain [m]



vessel_area = assemble(Constant(1)*(ds(1, domain = mesh) + ds(2, domain = mesh)))

N_a = len(artr)
N_v = len(vein)                                # Number of arteries and veins in box
N = N_a + N_v
vessel_length = N*ratio                        # [m] The total length of all vessels not including capillaries


print('Length of rectangle : ', ratio, ' m')


A_a = 2*plt.pi*R_a*N_a*1e-3                # Arterial surface area per meter cube [m]
A_v = 2*plt.pi*R_v*N_v*1e-3                # Venous surface area per meter        [m]
A = A_a + A_v                              # Total vessel area per meter          [m^2]
                                 

Q = ratio*1e3*flow_cube*1e-9               # total flow in m^3/s (from m*mm^2/s)
A_tot = A*ratio                            # total vessel suface area [m^2]


integrated_velocity = assemble(sqrt(inner(K/mu*grad(p_), K/mu*grad(p_)))*dx)  # mm/s*mm = mm^3/s
average_velocity = integrated_velocity/(A*1e6)*1e-3                           # [m/s]

print('Average Velocity = ', average_velocity, ' m/s')

print('Total flow = ', Q*1e6*60, ' mL/min')
print('Vessel area = ', A_tot, ' m^2')

R = Q**-1*(float(p_a) - float(p_v))        # Resistance in Pa*s/m^3
p_convert = 1e6*133*60.                    # converting from mmHg/(ml/min) to Pa *s /m^3

print('Total Resistance = ', R/p_convert, 'mmHg/(mL/min)')                         # Resistance in mmHg/(mL/min)
print('Total vessel length = ', vessel_length, ' m')


asgari_convert = 133**-1*1e12*60**-1



#print('IEG R = ', AU_r, 'mmHg/(mL/min)')

d_a = 24*1e-9                                  # Cleft size [m] (20 nm)
d_v = 31*1e-9

T_EF = 1e-6                                    # Cleft channel thickness [m]
W_EF = 5e-6
ieg = 20e-9

A_Asgari = W_EF*ieg                  # Cleft area in m^2. A_asgari*N_asgari = S*0.003
N_Asgari = A_tot*0.003/(A_Asgari)

R_Asgari = 95*asgari_convert

R_Asgari_tot = (1/R_Asgari*N_Asgari)**-1
#print('No clefts in Asgari: ', N_Asgari, 'Total cleft resistance from Asgari: ', R_Asgari_tot)



n_clefts_a = 3#(A_a*0.003)/(d_a*N_a)        # number of clefts per artery
n_clefts_v = 3#(A_v*0.003)/(d_v*N_v)        # s/S = 0.003 where s is average cleft surface and S is average enfeet/vessel surface. 
                                              # In 2D cross section: s = n_clefts*d where d ~ 20 nm. S = A/N, where a is surface area of vessels in the small block. 
                                              # n_clefts = A/ratio*0.003/d
                                              # See Asgari (2015) citations to Mathiisen et al (2010).




L_a = ratio                             # Length of one artery [m]
L_v = ratio                             # Length of one vein   [m]


R_cleft_a = 12*mu*T_EF/(d_a**3*L_a)     # Resistance of one cleft along entire artery [Pa/(m*s)]
R_cleft_v = 12*mu*T_EF/(d_v**3*L_v)     # Resistance of one cleft along entire vein   [Pa/(m*s)]

#print(n_clefts_a, n_clefts_v, L_a, L_v, R_cleft_a/p_convert, R_cleft_v/p_convert)

R_clefts_a = ((1/R_cleft_a)*N_a*n_clefts_a)**-1
R_clefts_v = ((1/R_cleft_v)*N_v*n_clefts_v)**-1


print('R_clefts_A = ', R_clefts_a/p_convert, 'R_clefts_V = ', R_clefts_v/p_convert)


print('Number of arteries: ', N_a, '  Number of veins: ', N_v,'\n')

print('Lumped resistance clefts_a + ECS + clefts_v = ', (R_clefts_a + R + R_clefts_v)/p_convert)


#print('Asgari flow rate = ', N*n_clefts*ratio*2*10**-14*60, ' mL/min')

File('pressure.pvd')<< p_













