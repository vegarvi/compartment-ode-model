from pylab import *
import scipy.optimize as sc
import scipy.signal as sig
import scipy.interpolate as sip
import pandas as pd
import IPython

frisk = True
if frisk:
    A = pd.read_csv('UMU 4356 Data_from_the_experiment_frisk.txt',decimal =',',delimiter='\t',names=['ICP','Needle','Flode','Code'])


    skip = 180          # skip first 160 seconds of ICP measurement

    p_data = A['ICP'][skip:]      # kPa
    q_data = A['Flode'][skip:]    # microL/sec
    Needle = A['Needle'][skip:]
    C = A['Code'][skip:]
    TIME = array([1595,4310,4800])-skip

else:
    A = pd.read_csv('UMU 5616 Data_from_the_experiment_iNPH.txt', decimal = '.', delimiter = '\t', names = ['Scaled ICP','Scaled IP','Flow','Unscaled IP', 'Pump frequency',  'Level', 'Status', 'Error', 'Balancing ICP', 'Balancing IP', 'Calibration ICP', 'Calibration IP'],skiprows=1)
    skip = 400          # skip first seconds of ICP measurement

    p_data = A['Scaled ICP'][skip:]      # kPa
    q_data = A['Flow'][skip:]    # microL/sec
    Needle = A['Scaled IP'][skip:]
    C = A['Status'][skip:]
    TIME = array([1295,3205,4195])-skip

def smoothen_flow(t_data,q_data):
    N_data = len(q_data)
    Q_exp = [trapz(q_data[:i], t_data[:i]) for i in range(1,N_data)]
    b, a = sig.butter(8,0.01)
    Q_vol = sig.filtfilt(b,a,Q_exp)
    Q_flow = (Q_vol[1:]-Q_vol[:-1])/(t_data[1]-t_data[0])
    Q = Q_flow*1e-3*1e-6   # from mircoL/sec to mL/sec to m^3/sec
    return Q


def optimize_params(p_data,q_data,TIME):
    N_data = len(p_data)
    t_data = linspace(0, N_data-1, N_data)
    
    Q = smoothen_flow(t_data,q_data)






    #pulsations related to cardiac cycle [Pa]
    # (not used at the moment)
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


    #pressure in dural sinus?  [Pa]
    pss = 4*133.33


    r_AG = 1e11   # CSF to arachnoid granulations

    p0 = 133.33*2


    p_legend = ['p_CSF'] 
    # A,C and V denote arterial, capillary and venous blood. a,c and v denote the corresponding PVS. 

    E = 2e5

    # physiological CSF and blood volumes + infusion. 

    Q_prod = 800*1e-6/24./60./60.   # 500 mL/day   (cm^3/day)
    Q_inf = sip.interp1d(t_data[:-2],Q) #Q_inf = 0.5*1e-6/60.            # 0.5 mL/min


    T = TIME[2]
    T_inf = TIME[0]
    T_inf_end = TIME[1]
    N = 32001

    



    def solve(P):
        t = linspace(0,T,N)
        p = zeros(N)
        CSF_start = average(p_data[0:100])/133.*1000
        E = exp(P[1])
        r_AG = exp(P[0])
        p[0] = CSF_start*133.33
        dt = t[1]-t[0]  

        for i in range(N-1):
            Q_in = Q_prod + Q_inf(t[i])          # Must have Q_inf(t) = 0 when there is no infusion. 
            Q_out = 1./r_AG*(p[i]-pss)#array([1./r_AG*(p[i,0]-pss), 0, 0, 1./r_V_AG*(p[i,3]-pss), 0, 0, 0, 0])

            p[i+1] = p[i] + dt/C(p[i],E,p0)*(Q_in - Q_out)
            
        p*= 1./133.33


        plot(t,p)
        draw()
        pause(0.01)
        p_end = p
        return t,p

    def functional(P):
        t,p = solve(P)
        N = len(t)
        p_func = sip.interp1d(t_data,p_data/133.*1000)  
        func = sum((p_func(t)-p))**2/N
        return func
            
            


    x0 = [log(r_AG),log(E)]
    bnds = [(xi*0.95,xi*1.05) for xi in x0]
    res = sc.minimize(functional, x0,options={'disp':True, 'maxiter':30},bounds = bnds)
    t,p = solve(res.x)
    return res, t, p

if __name__=='__main__':
    res, t, p = optimize_params(p_data,q_data,TIME)
    ioff()
    t_data = linspace(0,len(p_data)-1,len(p_data)) # 1 Hz
    plot(t_data,p_data/133.*1000,t,p)
    legend(['p_CSF_data','p_CSF_model'])
    xlabel('Time [s]')
    ylabel('Pressure [mmHg]')
    figure()
    plot(t_data[1:-1],1e6*smoothen_flow(t_data,q_data))
    legend('Q_infusion')
    xlabel('Time [s]')
    ylabel('Flow [mL/s]')
    show()
    

'''
functional = sum((90-p[:,1])**2)/N + sum((p_func(t)-p[:,0]))**2/N

res.x for iNPH:

res.x = array([27.90487736, 30.44602357, 20.55637488, 28.55750157, 17.97012865,
       27.40539972, 27.40579409, 25.47258903, 25.35744171, 27.15985949,
       28.10215859, 24.89575465, 25.66746445, 15.81402391, 11.97521788])

res.x for healthy:

res.x = array([27.64565145, 30.19334763, 20.52933878, 28.27941185, 18.25952455,
       27.65146048, 27.65155299, 25.6913622 , 25.57459917, 27.39126077,
       27.87078207, 25.10866926, 25.12033259, 15.97811175, 12.10025419])
'''
