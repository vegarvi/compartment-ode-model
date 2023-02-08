from pylab import *
import scipy.signal as sc
import pandas as pd



A = pd.read_csv('UMU 4356 Data_from_the_experiment.txt',decimal =',',delimiter='\t',names=['ICP','Needle','Flode','Code'])

ICP = A['ICP']
F = A['Flode']
Needle = A['Needle']
C = A['Code']

N = len(ICP)
f = 1.
t = linspace(0,(N-1)/f,N)
plot(t,ICP*100)
plot(t,F)
#plot(t,Needle*100)
plot(t,10*C)
show()

