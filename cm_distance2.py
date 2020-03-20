from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

table = '/home/elismarlosch/Simulations/galmer_tables/table_orbit1'

def r_cm(table_path, snapshot, galaxy):
    table = Table.read(table_path, snapshot)
    galaxy_mask = table['GAL'] == galaxy
    particle_mask = table['P_TYPE'][galaxy_mask] != 2
    m = table['MASS'][galaxy_mask][particle_mask]
    x = table['X'][galaxy_mask][particle_mask]
    y = table['Y'][galaxy_mask][particle_mask]
    z = table['Z'][galaxy_mask][particle_mask]
    mass = np.sum(m)
    x_cm = np.sum(x*m)/mass
    y_cm = np.sum(y*m)/mass
    z_cm = np.sum(z*m)/mass
    return np.array([x_cm, y_cm, z_cm])

def distance(table_path, snapshot):
    r_cm1 = r_cm(table_path, snapshot, 1)
    r_cm2 = r_cm(table_path, snapshot, 2)
    distance = np.sqrt(np.sum(np.square(r_cm1 - r_cm2)))
    return distance

#import pdb; pdb.set_trace()

d = np.zeros(71)
for i in range(1, 72):
    d[i - 1] = distance(table, i)

# r_cm1, r_cm2 = np.array([])
# for i in range(1, 72):
#     r_cm1 = np.concatenate((r_cm1, r_cm(table, i, 1)))
#     r_cm2 = np.concatenate((r_cm2, r_cm(table, i, 2)))

time = np.arange(0, 71, 1) * 50.0 #Myr

from matplotlib import rc
import matplotlib.font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

ax.set_title('distance between galaxies CMs', fontsize=25)
ax.tick_params(axis='both', tickdir='in', length=10, labelsize=20, width=1.5)
ax.set_xlabel(r'$time (Myr)$', fontsize=20)
ax.set_ylabel(r'$distance (kpc)$', fontsize=20)
ax.plot(time, d, '-.', color=colors[3])

plt.show()
