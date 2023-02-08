from pylab import *

n = array([2**18, 2**16, 2**16])
r = array([100.16, 97.42, 97.42])*1e-6/2


ratio_pvs = 1.26
A_vessel = pi*r**2
A_pvs = n*A_vessel*ratio_pvs

A_tot = sum(A_pvs)

print A_tot

p_convert = 1e6*133*60.        # converting from mmHg/(ml/min) to Pa *s /m^3                         
r_a = 1.14*p_convert              # Faghih

dp0 = 0.14*133.33
dp1 = 0.67*133.33

Q_pvs = array([dp0,dp1])/r_a

v_pvs = Q_pvs/A_tot

print v_pvs

