from pylab import *

def Lp(r0, i): 
    return 20*pow(0.5, i/3.54)*r0  # 3.54 for veins

def rp(r0, i): 
    return pow(0.5, i/3.54)*r0


def Resistance(r, mu, l, k): 
    RR = (8*mu/(pi*r**4))*(l / ((k**-4 -1) - (k**-2-1)**2/(log(k**-1))))
    return RR 
       

def Resistance_vessel(r, mu, l, k):     # VESSEL RESISTANCE: rename to Resistance and run to get venous tree resistance. 
    # Q = dp/l *pi*r**4/(8*mu)
    # R = dp/Q
    mu *= 3                # Blood viscosity
    R = 8*mu*l/(pi*r**4)
    return R


def f(k):
    return (1 / ((k**-4 -1) - (k**-2-1)**2/(log(k**-1))))

def Coeff(r,l):
    return l*(8*mu/(pi*r**4))


        


mu = 0.6951*10**-3 # Pa s 
k = 0.940721

n = 2*(2**30 + 2**28 + 2**28)
print(n)
r10 = 10e-6
l10 = 20*r10

p_factor = 60*133.*1e6



N = 11

Res = []

print('Generation  Radius  Resistance')

for i in range(N):
    Res.append(Resistance(r10, mu, l10, k)/n)
    print('%e %e %e %e %e'%(i, f(k), r10, Coeff(r10,l10), Res[-1]))
    r10 *= (pow(0.5, 1/3.54))**-1
    l10 = 20*r10
    n *= 0.5

print ('Resistance in mmHg/mL/min: ', sum(Res)/p_factor)



