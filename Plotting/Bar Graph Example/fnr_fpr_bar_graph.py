import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

plt.rcParams.update({'font.size': 13})

labels = ['150', '300', '375']

# # 2.12 - 0.04
# pcohe = [0.00923, 0.00281, 0.00147]
# ncohe = [0.01074, 0.00044, 0.00000]
# pside = [0.01224, 0.00280, 0.00107]
# nside = [0.00781, 0.00273, 0.00000]
# prand = [0.01327, 0.00239, 0.00128]
# nrand = [0.01283, 0.00050, 0.00161]
# pkirk = [0.02562, 0.01687, 0.00625]
# nkirk = [0.01917, 0.00738, 0.00000]

# # 3.98 - 0.08
# pcohe = [0.04288, 0.01084, 0.00639]
# ncohe = [0.02806, 0.00606, 0.00122]
# pside = [0.05140, 0.01169, 0.00488]
# nside = [0.03528, 0.00577, 0.00126]
# prand = [0.04263, 0.01003, 0.00616]
# nrand = [0.02764, 0.00847, 0.00142]
# pkirk = [0.06843, 0.02813, 0.01284]
# nkirk = [0.03691, 0.01708, 0.00558]

# # 6.01 - 0.12
# pcohe = [0.08477, 0.02934, 0.01319]
# ncohe = [0.04680, 0.01358, 0.00670]
# pside = [0.08900, 0.03038, 0.01348]
# nside = [0.04419, 0.01826, 0.00496]
# prand = [0.08074, 0.02344, 0.01662]
# nrand = [0.05764, 0.01626, 0.00756]
# pkirk = [0.10792, 0.06745, 0.03652]
# nkirk = [0.06139, 0.02607, 0.01806]

# 8.86 - 0.20
pcohe = [0.15968, 0.06830, 0.03654]
ncohe = [0.08570, 0.02652, 0.01204]
pside = [0.16455, 0.07121, 0.03816]
nside = [0.08650, 0.02134, 0.01405]
prand = [0.15250, 0.06547, 0.04141]
nrand = [0.07109, 0.02481, 0.01369]
pkirk = [0.19534, 0.11454, 0.06882]
nkirk = [0.12438, 0.03148, 0.01722]

x = np.arange(len(labels))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots()
ax.grid(b=True,which='both',axis='y', color='k')
ax.minorticks_on()
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.tick_params(axis='x', which='minor', bottom=False)
ax.set_axisbelow(True)

nrects1 = ax.bar(x - 3*width/2, np.array(nrand), width, label='Random Balanced',color='r',edgecolor='white', hatch='\\\\', linewidth=1)
nrects2 = ax.bar(x - 1*width/2, np.array(nkirk), width, label='Partial Kirkman',color='k',edgecolor='white', hatch='//', linewidth=1)
nrects3 = ax.bar(x + 1*width/2, np.array(ncohe), width, label='$\psi$-Optimal Balanced',color='m',edgecolor='white', hatch='xx', linewidth=1)
nrects4 = ax.bar(x + 3*width/2, np.array(nside), width, label='$\psi, \phi$-Optimal Balanced',color='b',edgecolor='white', hatch='--', linewidth=1)

ax.bar(x - 3*width/2, np.array(nrand), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x - 1*width/2, np.array(nkirk), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x + 1*width/2, np.array(ncohe), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x + 3*width/2, np.array(nside), width, color='none', edgecolor='k', linewidth=1)

prects1 = ax.bar(x - 3*width/2, -np.array(prand), width,color='r',edgecolor='white', hatch='\\\\', linewidth=1)
prects2 = ax.bar(x - 1*width/2, -np.array(pkirk), width,color='k',edgecolor='white', hatch='//', linewidth=1)
prects3 = ax.bar(x + 1*width/2, -np.array(pcohe), width,color='m',edgecolor='white', hatch='xx', linewidth=1)
prects4 = ax.bar(x + 3*width/2, -np.array(pside), width,color='b',edgecolor='white', hatch='--', linewidth=1)

ax.bar(x - 3*width/2, -np.array(prand), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x - 1*width/2, -np.array(pkirk), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x + 1*width/2, -np.array(pcohe), width, color='none', edgecolor='k', linewidth=1)
ax.bar(x + 3*width/2, -np.array(pside), width, color='none', edgecolor='k', linewidth=1)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('FPR              FNR', fontsize=25)
ax.set_xlabel('Number of measurements m', fontsize=25)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=18)
ax.legend()

# ax.set_yticklabels(['0.04','0.02','0.00','0.02','0.04'], fontsize=18)
# ax.set_yticks([-0.04,-0.02,0.00,0.02,0.04])
# ax.set_ylim(-0.04,0.04)

# ax.set_yticklabels(['0.08','0.06','0.04','0.02','0.00','0.02','0.04','0.06','0.08'], fontsize=18)
# ax.set_yticks([-0.08,-0.06,-0.04,-0.02,0.00,0.02,0.04,0.06,0.08])
# ax.set_ylim(-0.08,0.08)

# ax.set_yticklabels(['0.12','0.10','0.08','0.06','0.04','0.02','0.00','0.02','0.04','0.06','0.08','0.10','0.12'], fontsize=18)
# ax.set_yticks([-0.12,-0.10,-0.08,-0.06,-0.04,-0.02,0.00,0.02,0.04,0.06,0.08,0.10,0.12])
# ax.set_ylim(-0.12,0.12)

ax.set_yticklabels(['0.20','0.16','0.12','0.08','0.04','0.00','0.04','0.08','0.12','0.16','0.20'], fontsize=18)
ax.set_yticks([-0.20,-0.16,-0.12,-0.08,-0.04,0.00,0.04,0.08,0.12,0.16,0.20])
ax.set_ylim(-0.20,0.20)

# ax.set_yticklabels([str(abs(round(x,2))) for x in ax.get_yticks()], fontsize=18)
# ax.set_ylim(ax.get_yticks()[0],ax.get_yticks()[-1])

ax.axhline(y=0, linewidth=1, color='k')

fig.tight_layout()

plt.show()