from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt
from gas_to_the_center import r_cm_bp
import scipy.constants as cons


table_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
class MyTable():
    def __init__(self, table_path, snapshot, galaxy):
        self.table_path = table_path
        self.snapshot   = snapshot
        self.galaxy     = galaxy
        self.table      = Table.read(self.table_path, self.snapshot)
        self.masses, self.centers = self.binning(25)
        self.coords     = self.coordinates()


    def select_particles(self):
        x = self.table['X']
        y = self.table['Y']
        z = self.table['Z']
        r = np.array([x,y,z])
        radius_for_cm = 10
        cm_bp = r_cm_bp(self.table_path, self.snapshot, self.galaxy)
        mask = np.sqrt(np.sum(np.square(np.transpose(r) - cm_bp),axis=1)) < radius_for_cm
        return mask


    def binning(self, bins):
        mask = self.select_particles()
        x = self.table['X'][mask]
        y = self.table['Y'][mask]
        z = self.table['Z'][mask]
        coords = np.transpose(np.array([x, y, z]))
        masses, edges = np.histogramdd(coords, bins = bins, weights = self.table['MASS'][mask])
        # edges = np.asarray(edges)
        dbin = (edges[0][1] - edges[0][0]) / 2
        centers = np.array([i[:-1] for i in edges]) + dbin
        return masses, centers

    def coordinates(self):
        # coords = np.meshgrid(centers[0], centers[1], centers[2])
        coords = np.meshgrid(*self.centers)
        return np.asarray(coords)


    def distances(self, point):
        coords = self.coordinates()
        import pdb;pdb.set_trace()
        r = np.sqrt(np.sum([np.square(coords[i] - point[i]) for i in range(3)], 0))
#        dist = np.linalg.norm(coords1-coords2)
        return r


    def potential_energy(self):
        g = cons.G
        u = np.zeros((25, 25, 25))
        for i, j in np.ndenumerate(u):
            d = self.distances(self.coords[i])
            u[i] = -g*np.sum(self.masses/d)
        return u
