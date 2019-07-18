import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

x1, y1, error1 = np.loadtxt("mean_radius.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
#x2, y2 = np.loadtxt("actual_pos_dist.txt", usecols=(0,4), delimiter='\t', unpack='true')
X, Y, Z = np.loadtxt("actual_pos_dist.txt", usecols=(1,2,3), delimiter='\t', unpack='true')

plt.errorbar(x1,y1,yerr=error1)
plt.plot([1,50,100], [5,5,5])
plt.xlabel('# blocks')
plt.ylabel(r'$distance (\cdot a_0)$')


fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(X, Y, Z, c=Z, marker='.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.view_init(10, 30)
plt.show()

#plt.plot(x2,y2)



#plt.xscale('log')
#plt.yscale('log')

plt.show()
