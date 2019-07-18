import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit

def f(x,a,b):
    return a * np.exp(-b*x)

step, ene = np.loadtxt("inst_energy_liquid.txt",usecols=(0,1), delimiter='\t', unpack='true')

tmax=300

c1=0
c2=0
autocorrelation = np.zeros(tmax)


for t1 in range(0,100000):
    c1 += ene[t1]**2
    c2 += ene[t1]

sigma_2 = c1/len(step) - (c2/len(step))**2
#print(c1/len(step))
#print((c2/len(step))**2)
#print(sigma_2)
#print(len(step))


for t2 in tqdm(range(0,tmax)):
    c3=0
    c4=0
    c5=0
    delta=len(step)-t2
    for tt in range(0,delta):
        c3 += ene[tt]*ene[tt+t2]
        c4 += ene[tt]
        c5 += ene[tt+t2]

    autocorrelation[t2] = (c3/delta - c4*c5/(delta**2))/sigma_2
    #print(t2)
    #print(autocorrelation[t2])
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

#plt.yscale('log')

plt.title('Energy Autocorrelation - Liquid')
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Autocorrelation')
plt.show()
