import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit

def f(x,a,b):
    return a * np.exp(-b*x)

step, ene = np.loadtxt("inst_energy_solid.txt",usecols=(0,1), delimiter='\t', unpack='true')

tmax=300

c1=0
c2=0
autocorrelation = np.zeros(tmax)


for t in range(0,100000):
    c1 += ene[t]**2
    c2 +=ene[t]

sigma_2 = 1/len(step)*c1 - (1/len(step)*c2)**2
#print(sigma_2)


for t in tqdm(range(0,tmax)):
    c3=0
    c4=0
    c5=0
    delta=len(step)-t
    for tt in range(0,delta):
        c3 += ene[tt]*ene[tt+t]
        c4 += ene[tt]
        c5 += ene[tt+t]


    autocorrelation[t] = (1/delta*c3 -1/(delta**2)*c4*c5)/sigma_2


    #print(t)
    #print(autocorrelation[t])
    #print('\t')

x = np.arange(tmax)
#print(autocorrelation)

plt.plot(x, autocorrelation)

p_opt, p_cov = curve_fit(f, x, autocorrelation) #, bounds=([0,0,-1],[2,3,+3]))
y_fit = f(x,p_opt[0],p_opt[1])
plt.plot(x,y_fit) # plotting fitted function

print("optimized parameters [a,b] =")
print(p_opt)
print("parameters uncertainty =")
print(np.sqrt(np.diagonal(p_cov)))

tc = 1/p_opt[1]
print("Correlation time = ", tc)

plt.title('Energy Autocorrelation - Solid')
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Autocorrelation')
plt.show()
