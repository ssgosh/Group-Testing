import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import numpy as np

fact = [1,5,25,125]
factors = [59,60,61,62]

# Just some example data (random)
data = np.zeros((375,1000))

for i in range(4):
	data[0:375,125*i:125+125*i] = 1

for i in range(4):
	base = 500 + 125*i
	data[0:3*factors[i],base:base+factors[i]] = 1
	data[3*factors[i]:375,base+factors[i]:base+125] = 1

rows,cols = data.shape

fig = plt.figure() 
  
ax = fig.add_subplot(111)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.imshow(data, interpolation='nearest', 
                 extent=[0.5, 0.5+cols, 0.5, 0.5+rows],
                 cmap='tab20c')

for i in range(4):
	ax.add_patch( Rectangle((0.5+125*i, 0.5), 
	                        125, 375, 
	                        fc ='none',  
	                        ec ='k', 
	                        lw = 2) )
	plt.text(63+125*i, 188, str(fact[i]), horizontalalignment='center',verticalalignment='center', fontsize=15)

for i in range(4):
	ax.add_patch( Rectangle((500.5+125*i, 0.5+3*(125-factors[i])), 
	                        factors[i], 3*factors[i], 
	                        fc ='none',  
	                        ec ='k', 
	                        lw = 2) )
	plt.text(500.5+125*i + factors[i]/2.0, 375.5-3*factors[i]/2.0, str(factors[i]), horizontalalignment='center',verticalalignment='center', fontsize=15)
	ax.add_patch( Rectangle((500.5+125*i+factors[i], 0.5), 
	                        125-factors[i], 3*(125-factors[i]), 
	                        fc ='none',  
	                        ec ='k', 
	                        lw = 2) )
	plt.text(500.5+125*i+factors[i]+(125-factors[i])/2.0, 0.5+3*factors[i]/2.0, str(125-factors[i]), horizontalalignment='center',verticalalignment='center', fontsize=15)
plt.show()
