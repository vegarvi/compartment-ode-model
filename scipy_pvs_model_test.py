from pylab import *
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from IPython import embed


# 1.5 mL/min for 200 sec gives 3 mmHg pressure increase



plots = False
table = False
table2 = False
flow_data = False


mL_convert = 1e-6/60.    # converting from mL/min to m^3/s
model = 0

infusion_rate = 1.5*mL_convert    # 1.5 mL/min = 1.5*cm^3/(60s) = 1.5*1e-6m^3/60s = 1.5*1e-6/60. m^3/s
Q_prod = 0.33*mL_convert
Q_pulse_amp = 60*mL_convert      # 1 mL/sec = 60 mL/min


rest = 500          # time to let model reach steady state before infusion
infusion = 200      # duration of infusion test
T = 4800            # end time of simulation

t_inflow, flow = loadtxt('Inflow_data.txt')
t_ICP_data, ICP_data = loadtxt('ICP_data.txt')
ICP_data *= 1000/133.33
flow = 1e-3*60*flow*mL_convert                  # microL/sec to mL/min then convert

infusion_rate = interp1d(t_inflow, flow)

def Q_in(t, infusion = 1800, rest = 3600):
    Q_tot = Q_prod + Q_pulse_amp*sin(2*pi*t)
    #if type(t) == float or type(t) == int:
    #    if t > rest and t < rest+infusion:
    #        Q_tot += infusion_rate(t)
    #else:
    #    mask = (t > rest)*(t < rest + infusion)
     #   Q_tot += infusion_rate(t)*mask
    Q_tot += infusion_rate(t)

    return Q_tot

p_convert = 1e6*133*60.        # converting from mmHg/(ml/min) to Pa *s /m^3
p_AG =8.4*133.33
r_AG = 9.6*p_convert


p_base = 11.0*133.33
E = 0.3/(1e-6)
pr = 9*133.33

C_base = 1./(E*(p_base-pr))


def C(p,E,pr):
    if p < p_base:
        return C_base
    else:
        return 1./(E*(p-pr))
           
def f(p,t):
    dp = C(p,E,pr)**-1*(Q_in(t,infusion,rest) + (p_AG - p)/r_AG)

    return dp

p0 = 11*133.33


N = 10*T+1
t = linspace(0,T,int(N))
ICP, infodict = odeint(f,p0,t, full_output=True)
print(infodict['message'])

p = array(ICP[:,0])




# ALSO PLOT NO ABSORPTION BY CAPILLARIES
Q = 1./r_AG*(p_AG-p)





plot(t,p/133.33, t_ICP_data, ICP_data)
legend(['Model', 'Data'])
show()
'''
figure()
plot(t,-Q)
plot(t,Q_in(t,infusion, rest))
legend(['Outflow AG', 'Total inflow'])
vol_out = [trapz(Q[:i],t[:i]) for i in range(1,len(t))]
vol_in = [trapz(Q_in(t[:i],infusion,rest), t[:i]) for i in range(1,len(t))]
figure()
# volume is m^3 convert to cm^3 (mL) with 1e6
plot(t[1:],-1e6*array(vol_out), t[1:], 1e6*array(vol_in))
legend(['Integrated outflow', 'Integrated inflow'   ])
show()
'''

