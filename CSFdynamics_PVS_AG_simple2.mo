model CSFdynamics_PVS_AG_simple2

Real ICP;//unknown intracranial pressure
Real ICPmmHg;
Real Q_AG;//unknown absorption rate of CSF through AG
//Real Q_PVSa;//unknown absorption rate of CSF via PVS
Real Q_PVS;//unknown absorption rate of CSF via PVS

Real Qinf;//infusion rate during infusion test
Real Qpulse;//volume changes due to cardiac cycle
Real Qareg;//volume changes related to autoregulation
Real Q_in;//sum of inflow and production to the CSF space

Real CVP; //central venous pressure
Real Pd;//dural venous pressure
Real PL;//hypothetical venous pressure/lymph pressure??????
Real P0;//6*133;//P0 is a parameter related to the pressure in the extradural venous system, this parameter shifts the pressure-volume curve and do not demand max compliance at IPC = 0. Here it is assumed to be a function of Pd.
Real R_AG;
Real E;//elastance

parameter Real pi=3.1415;
parameter Real R_out=7.98e10; //= 10 mmHg/(ml/min) = 10 mmHg*min/ml = 10*133*1e6*60 Pa*s/m^3
parameter Real R_PVS = 5.0e14;//10*9.1e9; //= 1.14 mmHg/(ml/min) = 1.14 mmHg*min/ml = 1.14*133*1e6*60 Pa*s/m^3
//parameter Real R_AG=1.8e12; //= 226 mmHg/(ml/min) = 226 mmHg*min/ml = 226*133*1e6*60 Pa*s/m^3 -from Grzybowski.

parameter Real PVI=2e-5;//= 20 ml = 20*1e-6 m^3;//pressure volume index
parameter Real Qf=5e-9;//= 0.3 ml/min = 0.3*1e-6/60 m^3/s;

//Computations
//initial equation solving for ICP
initial equation
Qf=(1/R_AG)*(ICP-Pd) + (1/R_PVS)*(ICP-PL);

equation
R_AG*Q_AG=(ICP-Pd);
R_PVS*Q_PVS=(ICP-PL);

der(ICP)=E*(ICP-P0)*(Q_in-Q_AG-Q_PVS);

//R_PVS=R_out*R_AG/(R_AG-R_out);//8.35e10 Pa*s/m^3 
R_AG=R_out*R_PVS/(R_PVS-R_out); 

Qinf = if time > 5*60 and time < 25*60 then 1.5*1e-6/60 else 0;
//Qinf=0;//1.5*1e-6/60;
Qareg = 0.2*pi*0.67*1e-6*cos(2*pi/40*time);//corresponds to a pulsatile sourceterm with amplitude 0.1 ml and period 40 sec
Qpulse = pi*4e-6*cos(2*pi*time);//corresponds to a pulsatile sourceterm with amplitude 1 ml and period 1 sec

Q_in = Qf + Qinf + Qareg + Qpulse;

CVP=8*133;//+3*133*sin(2*pi*(8/60)*time);//sin related to respiration assuming 8 breath per min
Pd = CVP;//assume supine position
P0 = Pd - 2*133;//
PL = CVP + 2*133;
E = 1/(0.4343*PVI);

ICPmmHg = ICP/133.;

annotation(
   uses(Physiolibrary(version = "2.3.2-beta")));

end CSFdynamics_PVS_AG_simple2;