import networkx as nx
#import community
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def connectpoints(x1,y1,x2,y2):
    plt.plot([x1,x2],[y1,y2],'k-')



#partition = community.best_partition(G)  # compute communities

#G = nx.karate_club_graph()  # load a default graph
a, b =np.loadtxt("BestPath_3.txt", dtype=int, usecols=(0,1), delimiter='\t', unpack='true')
x, y = np.loadtxt("Map.txt", usecols=(1,2), delimiter='\t', unpack='true')  # compute graph layout

plt.figure()
plt.scatter(x,y)


for i in range(len(a)):
    tmp1=a[i]
    tmp2=b[i]
    connectpoints(x[tmp1],y[tmp1],x[tmp2],y[tmp2])


#pos = nx.spring_layout(G)

#print(pos)
#nx.draw(G, nx.spring_layout(G))


#plt.figure(figsize=(8, 8))  # image is 8 x 8 inches
#plt.axis('off')
#nx.draw_networkx_nodes(G, pos, node_size=100)#, cmap=plt.cm.RdYlBu, node_color=list(partition.values()))
#nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='b')

#plt.show(G)

'''
step, ene = np.loadtxt("inst_energy_gas.txt",usecols=(0,1), delimiter='\t', unpack='true')
press = np.loadtxt("inst_pressure_gas.txt",usecols=(1), delimiter='\t', unpack='true')


plt.plot(step, ene, label='Instant Energy')
plt.plot(step, press, label='Instant Pressure')
#plt.plot(mcstep, inst_m, label='Temperature')


plt.legend(loc='upper right')
'''

#x, gofr = np.loadtxt("output.epot.0",usecols=(0,1), delimiter='\t', unpack='true')
#x, y = np.loadtxt("Map.txt",usecols=(1,2), unpack='true')

#plt.scatter(x,y)

#f, axarr = plt.subplots(1,2,sharey=True)

#axarr[0].errorbar(x1,y1,yerr=error1)
#axarr[0].set_title('1 step')
#axarr[0].plot([1,50,100], [14.975790778311286, 14.975790778311286, 14.975790778311286])
#axarr[0].set(xlabel='# blocks')
#axarr[0].set(ylabel=r'$C[S(0),0]$')

#axarr[1].errorbar(x2,y2,yerr=error2)
#axarr[1].plot([1,50,100], [14.975790778311286, 14.975790778311286, 14.975790778311286])
#axarr[1].set_title('100 step')
#axarr[1].set(xlabel='# blocks')



#plt.xscale('log')
#plt.yscale('log')

plt.show()
