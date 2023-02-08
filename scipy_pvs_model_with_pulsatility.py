from pylab import *
from scipy.integrate import odeint
from IPython import embed

plots = False
table = False
table2 = False
flow_data = False


mL_convert = 1e-6/60.    # converting from mL/min to m^3/s
model = 0

infusion_rate = 1.5*mL_convert    # 1.5 mL/min = 1.5*cm^3/(60s) = 1.5*1e-6m^3/60s = 1.5*1e-6/60. m^3/s
Q_prod = 0.33*mL_convert
Q_pulse_amp = 200*mL_convert      # 1 mL/sec = 60 mL/min


rest = 1
infusion = 1800

def Q_in(t, infusion = 1800, rest = 3600):
    Q_tot = Q_prod + Q_pulse_amp*sin(2*pi*t)
    if type(t) == float or type(t) == int:
        if t> rest and t < rest+infusion:
            Q_tot += infusion_rate
    else:
        mask = (t > rest)*(t < rest + infusion)
        Q_tot += infusion_rate*mask
    

    return Q_tot

p_convert = 1e6*133*60.        # converting from mmHg/(ml/min) to Pa *s /m^3

                          
r_a = 1.14*p_convert              # Faghih


r_v = 1.75e-3*p_convert           # Faghih
r_V = 1.69*p_convert


r_cap = 125.31*p_convert
r_CP = 67*p_convert                    # Wild Guess ?
r_AG = 10.81*p_convert


p_AG =8.4*133.33
p_lymph = p_AG                 # Wild guess
p_CP = 0*133.33                
p_cap = 20*133.33


# REFERENCE MODEL
if model == 0: # only AG    #ref model: r = 8.6
    r_AG = 8.6*p_convert
    r_cap *= 1e9
    r_CP  *= 1e9
    r_V *= 1e9
    r_v *= 1e9
# MODIFICATION 2
if model == 2:
    r_a *= 1e9
    r_CP *= 1e9
# MODIFICATION 3
if model == 3:
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
    p_lymph = p_CP


r_lymph = r_V


p_base = 11.0*133.33
E = 0.2/(1e-6)
pr = 9*133.33

C_base = 1./(E*(p_base-pr))


def C(p,E,pr):
    if p < p_base:
        return C_base
    else:
        return 1./(E*(p-pr))
        






    
def f(p,t):
    # MODIFICATION 8
    p_pvs = (p_cap/r_cap + (r_lymph**-1 + r_a**-1)*p)*(r_cap**-1 + r_lymph**-1 + r_a**-1)**-1
    if p > p_cap:           # NO ABSORPTION BY CAPILLARIES
        p_pvs = p
    dp = C(p,E,pr)**-1*(Q_in(t,infusion,rest) + (p_CP - p)/r_CP + (p_AG - p)/r_AG + (p_pvs - p)/r_a + (p_pvs - p)/r_lymph)


    return dp

p0 = 11*133.33

T = 1802
N = 20*T+1
t = linspace(0,T,int(N))
ICP, infodict = odeint(f,p0,t, full_output=True)
print(infodict['message'])

p_CSF = array(ICP[:,0])
p_pvs = (p_cap/r_cap + (r_lymph**-1 + r_a**-1)*p_CSF)*(r_cap**-1 + r_lymph**-1 + r_a**-1)**-1
p = array([p_CSF, p_pvs]).transpose()



# ALSO PLOT NO ABSORPTION BY CAPILLARIES
Q = array([1./r_CP*(p_CP-p[:,0]), 1./r_AG*(p_AG-p[:,0]), 1./r_a*(p[:,1]-p[:,0])*(p[:,1] < p_cap), 1./r_cap*(p_cap-p[:,1])*(p[:,1] < p_cap), 1./r_lymph*(p[:,1]-p[:,0])*(p[:,1] < p_cap)])
Q = Q.transpose()

Q*= 1/mL_convert
idx0 = where(t>rest)[0][0]
idx1 = where(t>rest+infusion - 1)[0][0]


t -= 60*60
t *= 1./60


tref,pref,Qref = loadtxt('PressureFlow_ref.txt')



if plots:
    subplot(2,1,1)
    plot(t,p_CSF/133.33, t, p_pvs/133.33)
    plot(tref,pref/133.33,'--k')
    xlabel('Time [min]')
    ylabel('Pressure [mmHg]')
    legend(['CSF','PVS','ref'])
    axis([-2,60,5,25])
    subplot(2,1,2)

leg = ['crib','AG','PVS','cap','V']

if table:
    print('%d & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) \\\ '%(model, min(p_CSF[idx0:idx1])/133.33, max(p_CSF[idx0:idx1])/133.33, Q[idx0,1], Q[idx1,1], Q[idx0,2], Q[idx1,2], Q[idx0,0], Q[idx1,0], Q[idx0,3], Q[idx1,3]))

elif table2:
    print('%d %.2f %.2f %.2f & %.2f & %.2f \\\ '%(model, 1/mL_convert*(Q_prod+Q_in(0)), max(p_CSF)/133.33, p_convert**-1*(max(p_CSF[idx0:idx1])-min(p_CSF[idx0:idx1]))/Q_in(0), Q[idx0,3], Q[idx1,3]))
else:
    print('Clinical R_out = ', p_convert**-1*(max(p_CSF[idx0:idx1])-min(p_CSF[idx0:idx1]))/infusion_rate)
    print('min p ', min(p_CSF[idx0:idx1])/133.33, ' max p = ', max(p_CSF[idx0:idx1])/133.33)
    print('min pvs ', min(p_pvs[idx0:idx1])/133.33, ' max pvs = ', max(p_pvs[idx0:idx1])/133.33)  

for i in range(5):
    if plots:
        plot(t,Q[:,i])
    if flow_data:
        print(leg[i], '%.2f (%.2f)'%(Q[idx0,i], Q[idx1,i]))


plot(tref,Qref,'--k')
leg.append('ref')

if plots:
    legend(leg)
    axis([0,60,-1.85,0.3])
    xlabel('Time [min]')
    ylabel('Flow [mL/min]')


plot(t,p_CSF/133.33)
figure()
plot(t,Q[:,1])
plot(t,Q_in(t))
vol = [trapz(Q[:i,1],t[:i]) for i in range(1,len(t))]
figure()
plot(t[1:],vol)
show()


