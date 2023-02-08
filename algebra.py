from sympy import *
import sys
import pylab as plt
from IPython import embed


p_convert = 1e6*133*60.         # converting from mmHg/(ml/min) to Pa *s /m^3
mL_convert = 1e-6/60.           # converting from mL/min to m^3/s

def algebraic(model):
    p0 = Symbol('p0')
    p1 = Symbol('p1')

    pAG = Symbol('pAG')
    pcr = Symbol('pcr')
    pca = Symbol('pca')
    pL = Symbol('pL')

    rAG = Symbol('rAG')
    rcr = Symbol('rcr')
    rca = Symbol('rca')
    rL = Symbol('rL')
    rPVS = Symbol('rPVS')

    Q_prod = 0.33



    symbols = [pAG, pcr, pca, pL, rAG, rcr, rca, rL, rPVS]       
    r_a = 1.14*p_convert              # Faghih
    r_ECS = 79.78*p_convert           # Karl-Erik and distance between artery and vein in a 260 mum x 0.07 radius disc


    r_v = 1.75e-3*p_convert           # Faghih
    r_c = 32.24*p_convert


    r_cap = 125.31*p_convert                  # already multiplied with p_convert
    r_crib = 67*p_convert                    # Wild Guess ?
    r_AG = 10.81*p_convert


    p_AG = 8.8*133.33
    #p_AG = 11.5*133.33
    p_lymph = p_AG                 # Wild guess
    p_crib = 2*133.33                
    p_cap = 25*133.33

    #p_cap = 12.2*133.33             # Karens suggestion



    # REFERENCE MODEL
    if model == 0: # only AG    #ref model: r = 8.6
        r_AG = 8.6*p_convert
        r_cap *= 1e9
        r_crib  *= 1e9
        r_ECS *= 1e9
        r_c *= 1e9
    # MODIFICATION 2  # only AG w r = 10.81
    if model == 2:
        r_a *= 1e9
        r_crib *= 1e9
    # MODIFICATION 3
    if model == 3:    # elimiation of PVS
        r_a *= 1e9
    # MODIFICATION 4
    if model == 4:
        r_AG *= 1e9
    # MODIFICATION 5
    if model == 5:
        r_AG *= 2
    # MODIFICATION 6 and 8:
    if model == 6 or model == 8:
        r_a = 1.43e-3*p_convert
        r_c = 0.9e-3*p_convert
    # MODIFICATION 7:
    if model == 7:
        p_lymph = p_crib

    # MODIFICATION 9:
    if model == 9:
        Q_cap = 0.16/1e6/60.


    r_lymph = ((1./r_c + 1./r_ECS)**-1+r_v)


    pressures = plt.array([p_AG, p_crib, p_cap, p_lymph])



    resistances = plt.array([r_AG, r_crib, r_cap, r_lymph, r_a])


    values = list(pressures)+list(resistances)

    Q = Symbol('Q')

    if model == 8:
        F1 = Q + 1/rAG*(pAG - p0) + 1/rcr*(pcr - p0) + 1/rPVS*(p1 - p0) + 1/rL*(p1 - p0)
        F2 = 1/rL*(p0 - p1) + 1/rca*(pca - p1) + 1/rPVS*(p0 - p1)

    elif model == 9:
        F1 = Q + 1/rAG*(pAG - p0) + 1/rcr*(pcr - p0) + 1/rPVS*(p1 - p0)
        F2 = Q_cap + 1/rL*(pL - p1) + 1/rPVS*(p0 - p1)


    else:
        F1 = Q + 1/rAG*(pAG - p0) + 1/rcr*(pcr - p0) + 1/rPVS*(p1 - p0)
        F2 = 1/rL*(pL - p1) + 1/rca*(pca - p1) + 1/rPVS*(p0 - p1)



    sol0 = solve(F2,p1)[0]
    sol1 = solve(F1,p0)[0]

    F0 = F1.subs(p1,sol0)
    F1 = F2.subs(p0,sol1)

    P0 = (solve(F0,p0))[0]
    P1 = (solve(F1,p1))[0]


    replacements = [(s,v) for (s,v) in zip(symbols, values)]
    CSF = P0.subs(replacements)
    PVS = P1.subs(replacements)
    if model == 8:
        flux = [1/rAG*(pAG - CSF), 1/rPVS*(PVS - CSF) + 1/rL*(PVS - CSF), 1/rcr*(pcr - PVS), 1/rca*(pca - CSF)]
    elif model == 9:
        flux = [1/rAG*(pAG - CSF), 1/rPVS*(PVS - CSF), 1/rcr*(pcr - CSF), 0*rca*pcr]
    else:    
        flux = [1/rAG*(pAG - CSF), 1/rPVS*(PVS - CSF), 1/rcr*(pcr - CSF), 1/rca*(pca - PVS)]

    flux0 = plt.zeros(len(flux))
    flux1 = plt.zeros(len(flux))
    for i in range(len(flux)):
        flux[i] = flux[i].subs(replacements)
        flux0[i] = flux[i].subs(Q,0.33*mL_convert)/mL_convert
        flux1[i] = flux[i].subs(Q,1.83*mL_convert)/mL_convert


    CSF *= 1./133
    PVS *= 1./133
    CSF = CSF.subs(Q,Q/1e6/60.)  # SCALING
    PVS = PVS.subs(Q,Q/1e6/60.)
    CSF = CSF.subs(Q, Q + 0.33)
    PVS = PVS.subs(Q, Q + 0.33)
    flux = [1/rAG*(pAG - CSF), 1/rPVS*(PVS - CSF), 1/rcr*(pcr - CSF), 1/rca*(pca - CSF)]
    flux = [flux[i].subs(replacements) for i in range(len(flux))]
    R0 = float(CSF.args[1].subs(Q,1.0))
    R1 = float(PVS.args[1].subs(Q,1.0))
    #print('ICP SAS %.2f' % (CSF.subs(Q,0)))
    #print('FLOW PVS SAS %.2f' %(flux0[1]))
    #print('R PVS SAS %.2f'%(r_a/p_convert))
    #print('PVS pressure %.2f' %(CSF.subs(Q,1.5) + r_a*flux1[1]/p_convert))
    print('%d & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) & %.2f & %.2f \\\ '%(model, CSF.subs(Q,0), CSF.subs(Q,1.5), flux0[0], flux1[0], flux0[1], flux1[1], flux0[2], flux1[2], flux0[3], flux1[3], R0, R1))

for model in range(10):
    algebraic(model)


