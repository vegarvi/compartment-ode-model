from numpy import *

H_IEGa = 24e-9
H_IEGv = 31e-9
r_IEGa = H_IEGa/2
r_IEGv = H_IEGv/2
mu = 7e-4

L_EF = 15.7e-6
W_EF = 5e-6
T_EF = 1e-6
r_pvsa = 15e-6
r_pvsv = 20e-6

L_pvsa = 20*r_pvsa
L_pvsv = 20*r_pvsv


# number of AU along one artery or venule:
n_AUa = L_pvsa/W_EF
n_AUv = L_pvsv/W_EF
# number of gaps in on vessel
n_IEGa = (2*pi*r_pvsa/L_EF)*n_AUa
n_IEGv = (2*pi*r_pvsv/L_EF)*n_AUv


Ra = 8*mu*T_EF/(pi*r_IEGa**4*n_IEGa)
Rv = 8*mu*T_EF/(pi*r_IEGv**4*n_IEGv)

print("Ra", Ra)
print("Rv", Rv)
print("na", n_IEGa)
print("nv", n_IEGv)

convert = 133*60/(1e-6)

n_trees = 1.5625e9# alternative 1.67 m2/(area of vessels and veins)

Ra_trees = Ra/(n_trees*convert)
Rv_trees = Rv/(n_trees*convert)

print("Ra trees", Ra_trees)
print("Rv trees", Rv_trees)

