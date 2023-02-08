from pylab import *

#pressure in dural sinus?  [Pa]
pss = 4*133.33

#pulsations related to cardiac cycle [Pa]
def pulsation(t):
	return pi*sin(2*pi*t)*1e-6

#Compliance [m3/Pa] / elastance [1/m3]
def C(p,E,p0):
    return 1./(E*(p-p0))

N_C = 8 # number of compartments

# resistance between different compartments
# KHS: do you have references for all these?

# p = {p_CSF, p_A, p_C, p_V, p_a, p_c, p_v, p_ECS} 
# A,C and V denote arterial, capillary and venous blood. a,c and v denote the corresponding PVS. 


r_CSF_a = 8e11 #CSF -> PVSa
#r_CSF_v = 1e13 #CSF -> PVSv
r_CSF_ECS = 1e13 #CSF -> interstitial space

r_A_C = 8e8 # Arteries -> capillaries 
r_A_a = 1.5e12 # arteries -> PVSa


r_C_V = 1e8 # Capillaries -> veins
r_C_c = 1.3e12 # Capillaries -> PVSc
#r_C_ECS = 1.3e12 # Capillaries -> ECS

r_V_v = 1.3e12 # Veins -> PVSv


r_a_c = 1.8e11
r_a_ECS = 1.6e11

r_c_v = 1e12
r_c_ECS = 1e12

r_v_ECS = 1e11


r_AG = 2e11   # CSF to arachnoid granulations
r_V_AG = 1e7  # Intracranial venous blood to saggital sinus blood/AG


R_CSF = array([0, 0, 0, 0, r_CSF_a, 0, 0, r_CSF_ECS])
R_A = array([0, r_A_C, 0, r_A_a, 0, 0, 0])
R_C = array([0, r_C_V, 0, r_C_c, 0, 0])
R_V = array([0, 0, 0, r_V_v, 0])
R_a = array([0, r_a_c, 0, r_a_ECS])
R_c = array([0, r_c_v, r_c_ECS])
R_v = array([0, r_v_ECS])
R_ECS = array([0])

R_i = [R_CSF, R_A, R_C, R_V, R_a, R_c, R_v, R_ECS] # upper half of resistance matrix

Tau = zeros((N_C,N_C))

# create symmetrix matrix Tau_i,j = 1./R[i,j] for R \neq 0. if R = 0 (no connection between networks) set Tau = 0. 
for i in range(N_C):
    for j in range(i,N_C):
        if R_i[i][j-i] != 0:
            Tau[i,j] = 1./R_i[i][j-i]
            Tau[j,i] = Tau[i,j]             




p_legend = ['p_CSF', 'p_A', 'p_C', 'p_V', 'p_a', 'p_c', 'p_v', 'p_ECS'] 
# A,C and V denote arterial, capillary and venous blood. a,c and v denote the corresponding PVS. 


# no clue about these values 
#KHS: I can check if I find values! At least I can find for ECS and vessel walls.

E = array([1e3 for i in range(N_C)])
E[0] = 2e6
E[4] = 2e7
E[5] = 2e7
E[6] = 2e7
E[7] = 2e7
p0 = 133.33*array([2 for i in range(N_C)])

# physiological CSF and blood volumes + infusion. 

Q_prod = 500*1e-6/24./60./60.   # 500 mL/day   (cm^3/day)
Q_inf = 0.5*1e-6/60.            # 0.5 mL/min
B_in = 750*1e-6/60.             # 750 mL/min


Q_in = array([Q_prod, B_in, 0, 0, 0, 0, 0, 0])
Q_out = zeros(N_C)


T = 3000
T_inf = 1000
T_inf_end = 2000
N = 6001
t = linspace(0,T,N)
p = zeros((N,N_C))
p[0,:] = array([10, 85, 18, 10, 10, 10, 10, 10])*133.33
dt = t[1]-t[0]  

for i in range(N-1):
    if t[i] > T_inf and t[i]<T_inf_end:
        Q_in[0] = Q_prod + Q_inf
    elif t[i] >= T_inf_end:
        Q_in[0] = Q_prod
    Q_out = array([1./r_AG*(p[i,0]-pss), 0, 0, 1./r_V_AG*(p[i,3]-pss), 0, 0, 0, 0])

    p[i+1,:] = p[i,:] + dt/C(p[i,:],E,p0)*(Q_in-Q_out - sum([Tau[:,j]*(p[i,:]-p[i,j]) for j in range(N_C)],axis=0))

p*= 1./133.33
for i in range(8):
    plot(t,p[:,i])
legend(p_legend)
print(p_legend)
print(p[-1,:])
print(max(p[:,0])-p[-1,0])
show()


