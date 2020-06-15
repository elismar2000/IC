import matplotlib.pyplot as plt
import numpy as np
from pygadgetreader import *

#==========================
#Snapshots 0 and 1
#==========================

snap0 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_0000'
snap1 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_0001'

#==========================
#Reading ids and masses of gas in snapshot0. No star particles present.
#==========================

pid_gas0 = readsnap(snap0, 'pid', 'gas')

mass_gas0 = readsnap(snap0, 'mass', 'gas')

#==========================
#Reading ids and masses of gas and stars in snapshot1
#==========================

pid_gas1 = readsnap(snap1, 'pid', 'gas')
pid_star1 = readsnap(snap1, 'pid', 'star')

mass_gas1 = readsnap(snap1, 'mass', 'gas')
mass_star1 = readsnap(snap1, 'mass', 'star')

import pdb; pdb.set_trace()
#==========================
#Plotting stuff
#==========================

x = np.arange(0, len(mass_star1), 1)

plt.plot(x, mass_star1, '-', color='blue')
plt.show()

#110 partículas de gás perderam toda sua massa
#3745 partículas de gás perderam metade de sua massa
#Surgiram 3965 partículas de estrelas
#As partículas de gás que perderam todo seu gás foram transformadas em duas estrelas, metade da massa pra cada
#As partículas de gás que perderam metade do gás, foram transformadas em uma única estrela (toda massa de gás
#perdida virou a massa da estrela)
#A massa de todas as estrelas recém formadas é igual
