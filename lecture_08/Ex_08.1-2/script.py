import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
import subprocess
from shutil import *
from glob import glob
from tqdm import tqdm

#sigma=np.array([0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69])
#mu=np.array([0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89])

sigma=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0])
mu=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0])

min = 1000

with open('input.dat', 'r') as file:
  data = file.readlines()

#nblk = data[3]

for i in sigma:
    for j in mu:
        #print("sigma = ", i)
        #print("mu = ", j)
        data[1]=str(i)+"\n"
        data[2]=str(j)+"\n"
        with open('input.dat', 'w') as file:
            file.writelines( data )

        cmd= "./Variational_MC.exe"
        value = subprocess.call(cmd, shell = True)

        ene = np.loadtxt('Energy_sigma'+str(i)+"_mu"+str(j)+".txt",usecols=(1), delimiter='\t', unpack='true')
        if ene[99] < min:
            min = ene[99]
            s = i
            m = j


print("The parameters that minimize the energy are:")
print("sigma = ", s)
print("mu = ", m)
print("E_min = ", min)

iblk, e, err = np.loadtxt('Energy_sigma'+str(s)+"_mu"+str(m)+".txt",usecols=(0,1,2), delimiter='\t', unpack='true')

plt.errorbar(iblk, e, yerr=err)
plt.show()
