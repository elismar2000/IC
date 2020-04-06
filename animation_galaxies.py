from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit9'
snapshot = 1

t = Table.read(t_path, snapshot)

mask = t['P_TYPE'] != 2

cmap = plt.get_cmap('Accent')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, ax = plt.subplots()
plot, = ax.plot(t['X'], t['Y'], ',', color=colors[3])

title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.8, 'pad':5},
                transform=ax.transAxes, ha="center")
title.set_text('0 Myr')

def animate(snapshot):
    t = Table.read(t_path, snapshot)
    xlist = []
    ylist = []
    xlist.append(t['X'])
    ylist.append(t['Y'])

    time = (snapshot * 50) - 50

    title.set_text(u"{} Myr".format(time))

    plot.set_xdata(xlist)
    plot.set_ydata(ylist)
    return plot,

animation = FuncAnimation(fig, func=animate, frames=np.arange(1, 71, 1), interval=150)
plt.show()

animation.save('orbit9.mp4')
