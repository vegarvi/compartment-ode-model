from pylab import *


mL_convert = 1e-6/60.    # converting from mL/min to m^3/s

Q_in = 1.5*mL_convert    # 1.5 mL/min = 1.5*cm^3/(60s) = 1.5*1e-6m^3/60s = 1.5*1e-6/60. m^3/s
Q_prod = 0.33*mL_convert
Q_tot = Q_prod           # Start at baseline pressure

dt = 1
T = 8000
N = int(T/dt+1)

t = linspace(0,T,N)

p = zeros((N,2))



p_convert = 1e6*133*60.        # converting from mmHg/(ml/min) to Pa *s /m^3

                          
r_a = 1.14*p_convert              # Faghih
r_ECS = 2e15                   # Karl-Erik and distance between artery and vein in a 250 mum x 250 mum block?



r_v = 1.75e-3*p_convert           # Faghih
r_c_wall = 1e12
r_c = 32.24*p_convert


r_cap = (r_a+r_c_wall)
r_CP = 67*p_convert                    # Wild Guess ?
r_lymph = (r_c+r_v)
r_AG = 10.81*p_convert






p_AG = 5.9*133.33
p_lymph = 2*133.33             # Wild guess
p_CP = 2*133.33                # Wild guess
p_cap = 25*133.33




#p[0] = (1+r_out[0]/r_out[1])**-1*(p_out[0] + r_out[0]/r_out[1]*p_out[1] + Q_prod*r_out[0])  # steady state solution

p[0,0] = 12*133.33
p[0,1] = 12*133.33


def C(p,E,p0):
    return 1./(E*(p-p0))

E = 0.2/(1e-6)
p0 = 5*133.33

Q = zeros((N,5))


for i in range(N-1):    

    Q[i,:] = array([(p_CP - p[i,0])/r_CP, (p_AG - p[i,0])/r_AG, (p[i,1] - p[i,0])/r_a, (p_cap - p[i,1])/r_cap, (p_lymph - p[i,1])/r_lymph])

    if t[i]>=4000:
        Q_tot = Q_prod + Q_in
    if t[i] >= 6000:
        Q_tot = Q_prod

    p[i+1,0] = p[i,0] + dt/C(p[i,0],E,p0)*(Q_tot + (p_CP - p[i,0])/r_CP + (p_AG - p[i,0])/r_AG + (p[i,1] - p[i,0])/r_a)
    p[i+1,1] = p[i,1] + dt/C(p[i,1],E*10,p0)*( (p[i,0] - p[i,1])/r_a + (p_cap - p[i,1])/r_cap + (p_lymph - p[i,1])/r_lymph)


p_CSF = p[:,0]
p_pvs = p[:,1]

print('Clinical R_out = ', p_convert**-1*(max(p_CSF) - p_AG)/(Q_in+Q_prod), ' or ', p_convert**-1*(min(p_CSF) - p_AG)/Q_prod)
print('min p ', min(p_CSF)/133.33, ' max p = ', max(p_CSF)/133.33)    
subplot(2,1,1)
plot(t,p_CSF/133.33, t, p_pvs/133.33)
xlabel('Time [s]')
ylabel('Pressure [mmHg]')
legend(['CSF','pvs'])
subplot(2,1,2)
axis([0,2500,0,40])
Q*= 1/mL_convert
for i in range(5):
    plot(t,Q[:,i])
legend(['CP','AG','PVS','cap','lymph'])
axis([0,2500,-1.7,0.06])
xlabel('Time [s]')
ylabel('Flow [mL/min]')
show()
