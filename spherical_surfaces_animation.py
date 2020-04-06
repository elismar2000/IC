from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

mins = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins/orbit9_mins.txt')
x01, y01 = mins[0][1][0], mins[0][1][1]

t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit9'
snapshot = 1

t = Table.read(t_path, snapshot)

#mask1 = t['GAL'] == 1
gas = t['P_TYPE'] == 0

cmap = plt.get_cmap('Accent')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, ax = plt.subplots()
ax.set_xlim(-51.1, -49.0)
ax.set_ylim(-1.5, 1.0)

title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
                transform=ax.transAxes, ha="center")
title.set_text('0 Myr')

lines = []
plot1 = ax.plot(t['X'], t['Y'], ',', color=colors[3])[0]
plot2 = ax.plot(t['X'][gas], t['Y'][gas], '.', color='orange')[0]
lines.append(plot1)
lines.append(plot2)

u = np.mgrid[0:2*np.pi:40j]
radius = [1.0, 0.5, 0.1, 0.01]
for radius in radius:
    x1 = radius * np.cos(u) + x01
    y1 = radius * np.sin(u) + y01
    plot3 = ax.plot(x1, y1, '-', color='red')[0]
    lines.append(plot3)

def animate(snapshot):
    t = Table.read(t_path, snapshot)
    xlist = []
    ylist = []
    xlist.append(t['X'])
    ylist.append(t['Y'])
    xlist.append(t['X'][gas])
    ylist.append(t['Y'][gas])

    time = (snapshot * 50) - 50

    if snapshot < 72:
        x01, y01 = mins[snapshot - 1][1][0], mins[snapshot - 1][1][1]
    if snapshot >= 72:
        x01, y01 = mins[snapshot - 1][3][0], mins[snapshot - 1][3][1]

    ax.set_xlim(x01 - 1.1, x01 + 1.1)
    ax.set_ylim(y01 - 1.1, y01 + 1.1)
    title.set_text(u"{} Myr".format(time))

    u = np.mgrid[0:2*np.pi:40j]
    radius = [1.0, 0.5, 0.1, 0.01]
    for radius in radius:
        x1 = radius * np.cos(u) + x01
        y1 = radius * np.sin(u) + y01
        xlist.append(x1)
        ylist.append(y1)

    for lnum, line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum])

    return lines

animation = FuncAnimation(fig, func=animate, frames=np.arange(1, 71, 1), interval=500)
plt.show()

animation.save('inflow_orbit9.mp4')
