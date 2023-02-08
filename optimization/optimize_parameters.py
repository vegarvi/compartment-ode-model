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
    P_A = 100.
    P_C = 25.

    r_CSF_a = 8e11 #CSF -> PVSa
    #r_CSF_v = 1e13 #CSF -> PVSv
    r_CSF_ECS = 1e13 #CSF -> interstitial space

    r_A_C = 1e9 # Arteries -> capillaries 
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


    r_AG = 1e11   # CSF to arachnoid granulations
    r_V_AG = 1e7  # Intracranial venous blood to saggital sinus blood/AG

    def create_tau(R_i):
        r_CSF_a, r_CSF_ECS, r_A_C, r_A_a, r_C_V, r_C_c, r_V_v, r_a_c, r_a_ECS, r_c_v, r_c_ECS, r_v_ECS = R_i
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
        return Tau



    p_legend = ['p_CSF', 'p_A', 'p_C', 'p_V', 'p_a', 'p_c', 'p_v', 'p_ECS'] 
    # A,C and V denote arterial, capillary and venous blood. a,c and v denote the corresponding PVS. 


    # no clue about these values 
    #KHS: I can check if I find values! At least I can find for ECS and vessel walls.

    E0 = [1e3 for i in range(N_C)]
    E0[0] = 2e5
    E0[4] = 2e7
    E0[5] = 2e7
    E0[6] = 2e7
    E0[7] = 2e7
    p0 = 133.33*array([2 for i in range(N_C)])

    E = array(E0)

    # physiological CSF and blood volumes + infusion. 

    Q_prod = 800*1e-6/24./60./60.   # 500 mL/day   (cm^3/day)
    Q_inf = sip.interp1d(t_data[:-2],Q) #Q_inf = 0.5*1e-6/60.            # 0.5 mL/min
    B_in = 750*1e-6/60.             # 750 mL/min


    Q_in = array([Q_prod, B_in, 0, 0, 0, 0, 0, 0])
    Q_out = zeros(N_C)


    T = TIME[2]
    T_inf = TIME[0]
    T_inf_end = TIME[1]
    N = 32001
    Rij = [r_CSF_a, r_CSF_ECS, r_A_C, r_A_a, r_C_V, r_C_c, r_V_v, r_a_c, r_a_ECS, r_c_v, r_c_ECS, r_v_ECS, r_AG, r_V_AG]
    

    #Rij[0] = 8e8
    #Rij[2] = 7.15e8
    #Rij[4] = 8.5e7
    #Rij[7] = 6.5e8#12
    #Rij[8] = 9.6e8#14
    #Rij[9] = 1.4e7
    #Rij[11] = 1.4e7
    #Rij[13] = 6.4e7
        
    R0 = [log(x) for x in Rij]

    print(R0)
    

    p_end = zeros((N,N_C))
    '''
    ion()
    subplot(3,1,1)
    plot(t_data,p_data/133.*1000)
    axis([0, t_data[-1], min(p_data)*0.5/133*1e3, max(p_data)*1.5/133*1e3])
    subplot(3,1,2)
    plot(t_data,zeros(N_data))
    axis([0, t_data[-1], 0, 150])
    subplot(3,1,3)
    plot(t_data,zeros(N_data))
    axis([0, t_data[-1], 0, 10])
    show()
    '''

    def solve(P):
        R_all = [exp(x) for x in P[:len(R0)]]
        R_i = R_all[:-2]

        E = array(E0)
        E[0] = exp(P[-1])
        Tau = create_tau(R_i)
        r_AG = R_all[-2]
        r_V_AG = R_all[-1]



        t = linspace(0,T,N)
        p = zeros((N,N_C))
        CSF_start = average(p_data[0:100])/133.*1000
        p[0,:] = array([CSF_start, P_A, P_C, pss/133.33, CSF_start, CSF_start, CSF_start, CSF_start])*133.33
        dt = t[1]-t[0]  

        for i in range(N-1):
            Q_in[0] = Q_prod + Q_inf(t[i])          # Must have Q_inf(t) = 0 when there is no infusion. 

            Q_out = array([1./r_AG*(p[i,0]-pss), 0, 0, 1./r_V_AG*(p[i,3]-pss), 0, 0, 0, 0])
            p[i+1,:] = p[i,:] + dt/C(p[i,:],E,p0)*(Q_in-Q_out - sum([Tau[:,j]*(p[i,:]-p[i,j]) for j in range(N_C)],axis=0))
            
        p*= 1./133.33

        subplot(3,1,1)
        plot(t,p[:,0])
        subplot(3,1,2)
        plot(t,p[:,1])
        subplot(3,1,3)
        plot(t,p[:,3])
        draw()
        pause(0.01)
        p_end = p
        return t,p

    def functional(P):
        t,p = solve(P)
        N = len(t)
        p_func = sip.interp1d(t_data,p_data/133.*1000)  
        func = sum((P_A-p[:,1])**2)/N + sum((p_func(t)-p[:,0]))**2/N
        return func
            
            

    x0 = R0+[log(2e5)]#+E0
    
    bnds = [(xi*0.9,xi*1.1) for xi in x0]
    res = sc.minimize(functional, x0,options={'disp':True, 'maxiter':30},bounds = bnds)
    t,p = solve(res.x)
    return res, t, p

if __name__=='__main__':
    res, t, p = optimize_params(p_data,q_data,TIME)
    ioff()
    t_data = linspace(0,len(p_data)-1,len(p_data)) # 1 Hz
    plot(t_data,p_data/133.*1000,t,p[:,0])
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
