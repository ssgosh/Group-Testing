import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.ticker import NullFormatter
import math

m = np.zeros((15,35))
mat = np.loadtxt("kirkman_15_35.txt", dtype="int", delimiter=",")

cmap = ListedColormap(['#ffffff', '#777777'])

# Display matrix
plt.matshow(m,cmap=cmap)

ax = plt.gca()

plt.axhline(y=-0.5, color='#bbbbbb', linestyle='-', linewidth=4)
plt.axhline(y=14.5, color='#bbbbbb', linestyle='-', linewidth=4)

plt.axvline(x=-0.5, color='#bbbbbb', linestyle='-', linewidth=4)
plt.axvline(x=4.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=9.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=14.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=19.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=24.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=29.5, color='#bbbbbb', linestyle='-', linewidth=2)
plt.axvline(x=34.5, color='#bbbbbb', linestyle='-', linewidth=4)

# Minor ticks
ax.set_xticks(np.arange(-.5, 35, 1), minor=True)
ax.set_yticks(np.arange(-.5, 15, 1), minor=True)

# Major ticks
ax.set_xticks(np.arange(-.5, 35, 5))
ax.set_yticks(np.arange(-.5, 15, 15))
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

plt.box(False)

# Gridlines based on minor ticks
ax.grid(which='minor', color='#dddddd', linestyle='-', linewidth=1, zorder=0)

plt.rcParams["font.family"] = "Times New Roman"
def annotate_dim(ax,xyfrom,xyto,textx,texty,text=None):
    if text is None:
        text = str(np.sqrt( (xyfrom[0]-xyto[0])**2 + (xyfrom[1]-xyto[1])**2 ))

    ax.annotate("",xyfrom,xyto,arrowprops=dict(arrowstyle='<|-|>,head_length=0.6,head_width=0.3',color='#444444',lw=1),annotation_clip=False)
    ax.text(textx,texty,text,ha='center',va='center',fontsize=23,color='#444444')


annotate_dim(plt.gca(),[-0.5,-1],[4.5,-1],2,-1.5,'5')
annotate_dim(plt.gca(),[-1,-0.5],[-1,14.5],-1.8,7,'15')
annotate_dim(plt.gca(),[-0.5,15],[34.5,15],17,15.8,'35')

# plt.show()
plt.savefig('kirkman_15x35.eps', format='eps', bbox_inches='tight')