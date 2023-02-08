from pylab import *


n = array([2**i for i in list(range(7))+list(range(5,-1,-1))])      
d = array([23.97, 19.17, 15.28, 12.08, 9.46, 7.32, 8, 11.51, 14.53, 17.79, 21.45, 25.70, 30.77])*1e-6
r = d/2.0
L = array([1267.6, 930.3, 543.6, 302.3, 161.2, 154.7, 243.9, 473.9, 272.3, 426.6, 632.5, 844.2, 936.3])*1e-6

v = array([8.2, 6.41, 5.05, 4.03, 3.29, 2.75, 2.3, 1.11, 1.4, 1.86, 2.56, 3.57, 4.97])*1e-3 #m/s

deltaP = array([6.93, 5.87, 4.02, 2.7, 1.82, 2.35, 2.62, 1.27, 0.61, 0.89, 1.31, 1.78, 2.01])*133.33

wall = array([4.84, 4.25, 3.81, 3.49, 3.27, 3.14, 0.309, 1.15, 1.45, 1.78, 2.15, 2.57, 3.08])*1e-6


start_branch, end_branch = 0,False

if end_branch: 
    n = n[start_branch:end_branch]
    d = d[start_branch:end_branch]
    r = r[start_branch:end_branch]
    L = L[start_branch:end_branch]
    v = v[start_branch:end_branch]
    deltaP = deltaP[start_branch:end_branch]
    wall = wall[start_branch:end_branch]


add_wall = False
faghih = False

ratio_pvs = linspace(1.26,0.8,len(d))

#ratio_pvs = 0.8



A_vessel = pi*r**2
A_pvs = A_vessel*ratio_pvs

if add_wall:
    r += wall

r_pvs = sqrt(A_pvs/pi + r**2)


if faghih:
    r_pvs = r + 0.1e-6

print('Vessel radii :', r)

#print r_pvs/r
#print(A_pvs/A_vessel)

n_capillaries = 1e11
n_trees = int(n_capillaries/max(n))

area = 2*pi*r*L*n
#n_trees = 15/sum(area)                           # Based on a surface area of 15 m^2. 

#n_trees = 3378082                           # Based on a total of 750 mL/min blood flow, and velocities and areas in Paynes network. 


 

print('n_trees: ', n_trees)


mu = 0.7e-3

k = r/r_pvs
#print 'k',k
#R_f = 8*mu/(pi*r**4)*(L/( (k**-4 - 1) - (k**-2 - 1)**2/log(k**-1)))
R_f = 8*mu*L/pi*( r_pvs**4 - r**4 - (r_pvs**2 - r**2)**2/(log(r_pvs/r)) )**-1



faghih_factor = 133*60/(1e-6)
R_f = R_f/(faghih_factor*n_trees)
#print R_f



R_inner = sum(R_f/n)#[sum(R_f[i]/n[i] for i in range(0,len(d)))]

print('Microcirculation area: ', sum(area)*n_trees)
print('Total capillary length: ', sum(L*n)*n_trees)

print('Total resistance R = ', R_inner)




'''
A = 3e-3*3e-3

mu = 7e-4

K = n*pi*r**4/(8*A)
print('K', K)

R_Q = deltaP/(v*n*pi*r**2)
print('R_Q', R_Q)

R_K = mu*L/(K*A)
print('R_K', R_K)

K_pvs = n*pi/(8*A)*(r_pvs**4-r**4 - (r_pvs**2-r**2)**2/log(r_pvs/r))
print('K_pvs', K_pvs)

R_Kpvs = mu*L/(K_pvs*A)
print('R_Kpvs', R_Kpvs)
 '''
