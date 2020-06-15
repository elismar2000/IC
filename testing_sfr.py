from __future__ import print_function

import numpy as np
from pygadgetreader import *

snap5 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_0005'
snap6 = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_0006'

ids5 = readsnap(snap5, 'pid', 'star')
ids6 = readsnap(snap6, 'pid', 'star')

new_born = np.isin(ids6, ids5)

mask = new_born == False

mass6 = readsnap(snap6, 'mass', 'star')

import pdb; pdb.set_trace()
