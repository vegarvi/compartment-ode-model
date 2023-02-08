from pylab import *

# Eq (3)
def Lp(r0, i): 
    return 20*pow(0.5, i/3.67)*r0 

# Eq (2)    
def rp(r0, i): 
    return pow(0.5, i/3.67)*r0 

# Eq (4)
def Resistance(r, mu, l, k): 
    RR = (8*mu/(pi*r**4))*(l / ((k**-4 -1) - (k**-2-1)**2/(log(k**-1))))
    return RR 


def f(k):
    return (1 / ((k**-4 -1) - (k**-2-1)**2/(log(k**-1))))

def Coeff(r,l):
    return l*(8*mu/(pi*r**4))

# For further studies
# VESSEL RESISTANCE: rename to Resistance and run to get arterial tree resistance. 
def Resistance(r, mu, l, k):    
    mu *= 3                               # Blood viscosity
    R = 8*mu*l/(pi*r**4)
    return R


rho = 993.0      # kg/m^3 # density
nu = 7*1e-7      # kinematic viscosity
mu = rho*nu      # Pa s CSF/WATER dynamic VISCOSITY
k = 0.6652 

r0_MCA = 1.5e-3         # Table 1
r0_ACA = 1e-3
r0_PCA2 =1e-3

l_MCA = 51e-3           # Table 1
l_ACA = 45e-3
l_PCA2 = 60e-3

p_factor = 60*133.*1e6   # To convert resistance to mmHg/mL/min. 


# LISTS OF TOTAL (PARALELL) RESISTANCE OF A GIVEN GENERATION IN THE TREE
# GEN 16 - 28 FOR MCA. GEN 18 - 30 FOR ACA AND PCA2
Res_MCA = []
Res_ACA = []
Res_PCA2 = []

N = 30


print('Gen', ' r_MCA', '     r_ACA', '    r_PCA2')

for i in range(16,N):
    r_MCA = rp(r0_MCA,i)
    r_ACA = rp(r0_ACA,i)
    r_PCA2 = rp(r0_PCA2,i)
              
    l_MCA = Lp(r0_MCA,i)    
    l_ACA = Lp(r0_ACA,i)
    l_PCA2 = Lp(r0_PCA2,i)
    


    # ADD TO MCA ONLY FOR 18 - 29:
    if i>= 18:
        Res_MCA.append(Resistance(r_MCA, mu, l_MCA, k)/2**i)

    # ADD TO ACA, PCA2 ONLY FOR 16 - 27:
    if i < N-2:                                                
        Res_ACA.append(Resistance(r_ACA, mu, l_ACA, k)/2**i)
        Res_PCA2.append(Resistance(r_PCA2, mu, l_PCA2, k)/2**i)
    
    
    '''
    print( '%.2d %10.4e %10.4e %10.4e %10.4e %10.4e'%(i, f(k), 2*r_MCA, Coeff(r_MCA, l_MCA), Resistance(r_MCA, mu, l_MCA, k), Resistance(r_MCA, mu, l_MCA, k)/2**i))     
    '''      
    # PRINT OUT A NICE TABLE OF RADII
    if i < 18:
        print( '%.2d %10s %10.4e %10.4e'%(i, 'N/A', r_ACA, r_PCA2)) 

    elif 18 <= i < N-2:
        print( '%.2d %10.4e %10.4e %10.4e'%(i, r_MCA, r_ACA, r_PCA2))       
    else:  
        print( '%.2d %10.4e %10s %10s'%(i, r_MCA, 'N/A', 'N/A'))   

    print('\n PCA2 Resistance: ', Resistance(r_PCA2, mu, l_PCA2, k), Resistance(r_PCA2, mu, l_PCA2, k)/2**i)
    print('r = %g    Coeff = %g,   f(k) = %g '%( r_PCA2, Coeff(r_PCA2, l_PCA2), f(k)))
    print()
    


# Radii and lengths of final generations 30, 28, 28:
r_MCA = rp(r0_MCA,30)
r_ACA = rp(r0_ACA,28)
r_PCA2 = rp(r0_PCA2,28)

l_MCA = Lp(r0_MCA,30)
l_ACA = Lp(r0_ACA,28)
l_PCA2 = Lp(r0_PCA2,28)

# print RESISTANCE OF FINAL GAPS OF PRECAPILLARIES

print(r_MCA, Resistance(r_MCA, mu, l_MCA, r_MCA/(r_MCA+100e-9))/2**30)
print(r_ACA, Resistance(r_ACA, mu, l_ACA, r_ACA/(r_ACA+100e-9))/2**28)
print(r_PCA2, Resistance(r_PCA2, mu, l_PCA2, r_PCA2/(r_PCA2+100e-9))/2**28)


Res_MCA.append(Resistance(r_MCA, mu, l_MCA, r_MCA/(r_MCA+100e-9))/2**30)
Res_ACA.append(Resistance(r_ACA, mu, l_ACA, r_ACA/(r_ACA+100e-9))/2**28)
Res_PCA2.append(Resistance(r_PCA2, mu, l_PCA2, r_PCA2/(r_PCA2+100e-9))/2**28)
print( 'N  %10.4e %10.4e %10.4e'%(r_MCA, r_ACA, r_PCA2))

#print( '%.2d %10.4e %10.4e %10.4e %10.4e %10.4e'%(i, r_MCA/(r_MCA+100e-9), 2*r_MCA, Coeff(r_MCA, l_MCA), Resistance(r_MCA, mu, l_MCA, k), Resistance(r_MCA, mu, l_MCA, k)/2**i))  

# SUM OVER ALL GENERATIONS IN A GIVEN TREE. PARALELL FLOW BETWEEN THE TREE TREES
# TWO SYMMETRIC TREES 

R_MCA = sum(Res_MCA)*0.5
R_ACA = sum(Res_ACA)*0.5
R_PCA2 = sum(Res_PCA2)*0.5


#Res = 0.5*(1./sum(Res_MCA) + 1./sum(Res_ACA) + 1./sum(Res_PCA2))**-1
Res = (1/R_MCA + 1/R_ACA + 1/R_PCA2)**-1
print('%e %e %e' %(R_PCA2, R_MCA, R_ACA))

#CONVERTING TO MMHG/ML/MIN GIVES:
print ('\nResistance in mmHg/mL/min: ', Res/p_factor)
    





