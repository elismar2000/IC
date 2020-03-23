from astropy.table import Table
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

table = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit1'

t = Table.read(table, 7)

x01, y01 = -12.9484, -7.4092
x02, y02 = 12.4350, 7.3115

mask1 = t['GAL'] == 1
mask2 = t['GAL'] == 2

cmap = plt.get_cmap('Accent')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(t['X'][mask1], t['Y'][mask1], t['Z'][mask1], ',', zdir='z', color=colors[3])
# ax.plot(t['X'][mask2], t['Y'][mask2], t['Z'][mask2], ',', zdir='z', color='orange')

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t['X'][mask1], t['Y'][mask1], ',', color=colors[3])
ax.plot(t['X'][mask2], t['Y'][mask2], ',', color=colors[1])

u = np.mgrid[0:2*np.pi:40j]
radius = [1.0, 0.5, 0.1, 0.01]
for radius in radius:
    x1 = radius * np.cos(u) + x01
    y1 = radius * np.sin(u) + y01
    x2 = radius * np.cos(u) + x02
    y2 = radius * np.sin(u) + y02
    ax.plot(x1, y1, '-', color='red')
    ax.plot(x2, y2, '-', color='red')

plt.show()

#AINDA FALTA:
#(a) Marcar as partículas de gás
#(b) Fazer isso para vários snapshots
#(c) Criar uma animação
