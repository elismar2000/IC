import matplotlib.pyplot as plt
import numpy as np
import pickle
from tables import *

orbit1_mins = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins_GalMer/orbit1_mins.txt'
mins_galmer_like_sim = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins_Gadget/galmer-like_sim_test2_minima.h5'


def distance(r1, r2):
    dist = np.sqrt(np.sum(np.square(r1[i] - r2[i]) for i in range(3)))
    return dist


with open(orbit1_mins, 'rb') as f:
    table = pickle.load(f)
    galmer_two_mins = np.sum(np.array([(type(table[i, 1]) == np.ndarray) for i in range(len(table[:, 1]))]))
    galmer_simulation_distance = np.array([distance(table[i, 1], table[i, 2]) for i in range(galmer_two_mins)])


gadget_simulation_h5file = open_file(mins_galmer_like_sim, 'a')
gadget_simulation_table = gadget_simulation_h5file.root.potential.readout

xmin1 = np.array([j['xmin1'] for j in gadget_simulation_table.iterrows()])
ymin1 = np.array([j['ymin1'] for j in gadget_simulation_table.iterrows()])
zmin1 = np.array([j['zmin1'] for j in gadget_simulation_table.iterrows()])
coords_min1 = np.vstack((xmin1, ymin1, zmin1))

xmin2 = np.array([j['xmin2'] for j in gadget_simulation_table.iterrows()])
ymin2 = np.array([j['ymin2'] for j in gadget_simulation_table.iterrows()])
zmin2 = np.array([j['zmin2'] for j in gadget_simulation_table.iterrows()])
coords_min2 = np.vstack((xmin2, ymin2, zmin2))

gadget_two_mins = np.sum(np.isnan(xmin1) == False)
gadget_simulation_distance = np.array([distance(coords_min1[:, i], coords_min2[:, i]) for i in range(gadget_two_mins)])

########################################################
time_galmer = np.arange(0, galmer_two_mins * 50, 50)
time_gadget = np.arange(0, gadget_two_mins * 50 * 0.98, 50)

fig = plt.figure()
ax = fig.add_subplot(111)

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

ax.tick_params(axis='both', tickdir='in', length=10, labelsize=15, width=1.5)
ax.set_xlabel(r'$time (Myr)$', fontsize=15)
ax.set_ylabel(r'$distance\ between\ galaxies\ (kpc)$', fontsize=15)
ax.plot(time_galmer, galmer_simulation_distance, '-.', color=colors[3], label='GalMer simulation')
ax.plot(time_gadget, gadget_simulation_distance, '-.', color=colors[5], label='Equivalent simulation but using Gadget-3')

plt.legend()
plt.show()
