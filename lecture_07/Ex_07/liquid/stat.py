import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit

def error(sum, sum2, iblk):
    return np.sqrt((sum2/iblk - (sum/iblk)**2)/iblk)

step, ene = np.loadtxt("inst_energy_liquid.txt",usecols=(0,1), delimiter='\t', unpack='true')
press = np.loadtxt("inst_pressure_liquid.txt",usecols=(1), delimiter='\t', unpack='true')


nsteps = 100000
iu=0
ip=1
glob_av = np.zeros(2)
glob_av2 = np.zeros(2)
blk_av = np.zeros(2)
error_u = np.zeros(100)
error_p = np.zeros(100)
x = np.arange(50,5001,50)
x[0]=10


for i in tqdm(range(len(x))):
    N_blocks = nsteps//x[i]
    glob_av[iu]=0
    glob_av[ip]=0
    glob_av2[iu]=0
    glob_av2[ip]=0
    #glob_av.fill(0)
    #glob_av2.fill(0)

    for j in range(N_blocks):
        blk_av[iu] = 0
        blk_av[ip] = 0

        for k in range(x[i]):
            #print(ene[block_size*j+k])
            blk_av[iu] += ene[x[i]*j+k]
            blk_av[ip] += press[x[i]*j+k]

        glob_av[iu] += blk_av[iu]/x[i]
        glob_av[ip] += blk_av[ip]/x[i]
        glob_av2[iu] += (blk_av[iu]/x[i])**2
        glob_av2[ip] += (blk_av[ip]/x[i])**2

    error_u[i] = error(glob_av[iu], glob_av2[iu], N_blocks)
    #print(glob_av[iu]**2)
    #print(glob_av2[iu])
    #print(error_u[i])
    error_p[i] = error(glob_av[ip], glob_av2[ip], N_blocks)


#plt.figure()
plt.plot(x,error_u)
#plt.figure()
plt.plot(x,error_p)

plt.title('Statistical Uncertainties')
plt.xlabel('Block Size')
plt.ylabel('Error')
plt.show()
