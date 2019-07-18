import matplotlib
import matplotlib.pyplot as plt
import numpy as np


ekin, err1 = np.loadtxt("ave_ekin.out", usecols=(0,1), delimiter=' ', unpack='true')
epot, err2 = np.loadtxt("ave_epot.out", usecols=(0,1), delimiter=' ', unpack='true')
etot, err3 = np.loadtxt("ave_etot.out", usecols=(0,1), delimiter=' ', unpack='true')
pr, err4 = np.loadtxt("ave_pr.out", usecols=(0,1), delimiter=' ', unpack='true')


plt.errorbar(ekin, yerr=err1, label='Kinetic Energy')
plt.errorbar(epot, yerr=err2, label='Potential Energy')
plt.errorbar(etot, yerr=err3, label='Total Energy')
plt.errorbar(pr, yerr=err4, label='Pressure')

plt.legend(loc='upper right')
