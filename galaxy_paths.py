import pickle
import numpy as np
import matplotlib.pyplot as plt

with open('orbit1_mins.txt', 'rb') as f:
    table1 = pickle.load(f)

x1_1 = np.zeros(71)
y1_1 = np.zeros(71)
z1_1 = np.zeros(71)
x2_1 = np.zeros(71)
y2_1 = np.zeros(71)
z2_1 = np.zeros(71)

for i in range(17):
    x1_1[i] = table1[i, 1][0]
    y1_1[i] = table1[i, 1][1]
    z1_1[i] = table1[i, 1][2]
    x2_1[i] = table1[i, 2][0]
    y2_1[i] = table1[i, 2][1]
    z2_1[i] = table1[i, 2][2]

r1_1 = np.array([x1_1, y1_1, z1_1])
r2_1 = np.array([x2_1, y2_1, z2_1])
###########################################################

with open('orbit5_mins.txt', 'rb') as f:
    table5 = pickle.load(f)

x1_5 = np.zeros(71)
y1_5 = np.zeros(71)
z1_5 = np.zeros(71)
x2_5 = np.zeros(71)
y2_5 = np.zeros(71)
z2_5 = np.zeros(71)

for i in range(29):
    x1_5[i] = table5[i, 1][0]
    y1_5[i] = table5[i, 1][1]
    z1_5[i] = table5[i, 1][2]
    x2_5[i] = table5[i, 2][0]
    y2_5[i] = table5[i, 2][1]
    z2_5[i] = table5[i, 2][2]

r1_5 = np.array([x1_5, y1_5, z1_5])
r2_5 = np.array([x2_5, y2_5, z2_5])
###########################################################

with open('orbit9_mins.txt', 'rb') as f:
    table9 = pickle.load(f)

x1_9 = np.zeros(71)
y1_9 = np.zeros(71)
z1_9 = np.zeros(71)
x2_9 = np.zeros(71)
y2_9 = np.zeros(71)
z2_9 = np.zeros(71)

for i in range(71):
    x1_9[i] = table9[i, 1][0]
    y1_9[i] = table9[i, 1][1]
    z1_9[i] = table9[i, 1][2]
    x2_9[i] = table9[i, 2][0]
    y2_9[i] = table9[i, 2][1]
    z2_9[i] = table9[i, 2][2]

r1_9 = np.array([x1_9, y1_9, z1_9])
r2_9 = np.array([x2_9, y2_9, z2_9])
#############################################################

which_plot = 1

cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

from matplotlib import rc
import matplotlib.font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('text', usetex=True)

#Pra plotar a distância entre as galáxias
if which_plot == 1:
    dist1 = np.array([np.linalg.norm(r1_1[:, i] - r2_1[:, i]) for i in range(71)])
    dist5 = np.array([np.linalg.norm(r1_5[:, i] - r2_5[:, i]) for i in range(71)])
    dist9 = np.array([np.linalg.norm(r1_9[:, i] - r2_9[:, i]) for i in range(71)])

    time = np.arange(0, 71, 1)*50 #Myr
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    ax.plot(time, dist1, '-.', color=colors[3], label='Orbit type 1', linewidth=2.0)
    ax.plot(time, dist5, '-.', color=colors[1], label='Orbit type 5', linewidth=2.0)
    ax.plot(time, dist9, '-.', color=colors[5], label='Orbit type 9', linewidth=2.0)
    ax.set_xlabel(r'$Myr$', fontsize=30)
    ax.set_ylabel(r'$kpc$', fontsize=30)
    ax.axvline(350, color='k', linestyle='dashed', linewidth=1.0)
    ax.axvline(400, color='k', linestyle='dashed', linewidth=1.5)
    ax.axvline(450, color='k', linestyle='dashed', linewidth=2.0)
    ax.axvline(1250, color='k', linestyle='dashed', linewidth=1.5)
    ax.tick_params(axis='both', tickdir='in', length=10, labelsize=20, width=1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)

    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    axins = zoomed_inset_axes(ax, 2.5, loc=10)
    axins.plot(time, dist1, '-.', color=colors[3], label='Orbit type 1')
    axins.plot(time, dist5, '-.', color=colors[1], label='Orbit type 5')
    axins.plot(time, dist9, '-.', color=colors[5], label='Orbit type 9')
    axins.tick_params(axis='both', tickdir='in', length=4, labelsize=20, width=1.0)
    plt.yticks(visible=False)
    plt.xticks(visible=False)

    from mpl_toolkits.axes_grid1.inset_locator import mark_inset
    mark_inset(ax, axins, loc1=2, loc2=4, fc='none', ec='0.5')

    x1, x2, y1, y2 = 280, 500, 10, 35
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    #plt.tight_layout()
    ax.legend(loc=1, fontsize=20, fancybox=True, framealpha=0.9, shadow=True)
    plt.show()

#Pra plotar a trajetória das galáxias
if which_plot == 2:
    fig, axs = plt.subplots(1, 3, figsize=(20, 20))

    axs[0].plot(x1_1, z1_1, '-')
    axs[0].plot(x2_1, z2_1, '-')
    axs[0].set_title('orbit1')
    axs[1].plot(x1_5, z1_5, '-')
    axs[1].plot(x2_5, z2_5, '-')
    axs[1].set_title('orbit5')
    axs[2].plot(x1_9, z1_9, '-')
    axs[2].plot(x2_9, z2_9, '-')
    axs[2].set_title('orbit9')
    plt.show()

if which_plot == 3:
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(x1_9, y1_9, z1_9, zdir='z')
    ax.plot(x2_9, y2_9, z2_9, zdir='z')
    ax.set_aspect('equal')
    plt.show()
