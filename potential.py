from typing import List

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cons
from astropy.table import Table
from scipy.ndimage import gaussian_filter

from gas_to_the_center import r_cm_bp
from particle_positions import position, position_all_particles


class MyTable:
    def __init__(self, table_path: str, snapshot: int, galaxy: int, n_bins: int = 25):
        self.table_path = table_path
        self.snapshot = snapshot
        self.galaxy = galaxy
        self.table = Table.read(self.table_path, self.snapshot)
        self.masses, self.centers, self.d_bin = self.binning(n_bins)
        self.coords = self.coordinates()
        self.n_bins = n_bins

    def select_particles(self):
        x = self.table['X']
        y = self.table['Y']
        z = self.table['Z']
        r = np.array([x, y, z])
        radius_for_cm = 20
        cm_bp = r_cm_bp(self.table_path, self.snapshot, self.galaxy)
        mask = np.sqrt(np.sum(np.square(np.transpose(r) - cm_bp), axis=1)) < radius_for_cm
        return mask


    def select_particles_cube(self, coords):
#        width = (self.coords[0, 1, 1, 1] - self.coords[0, 1, 0, 1]) / 2
        width = self.d_bin * 1.25
        mask = (self.table['X'] >= (coords[0] - width)) & (self.table['X'] <= (coords[0] + width))\
        & (self.table['Y'] >= (coords[1] - width)) & (self.table['Y'] <= (coords[1] + width))\
        & (self.table['Z'] >= (coords[2] - width)) & (self.table['Z'] <= (coords[2] + width))
#        import pdb; pdb.set_trace()
#        print('width = ', width)
        return mask


    def binning(self, bins: int, mask = None):
        if mask is None:
            x = self.table['X']
            y = self.table['Y']
            z = self.table['Z']
            w = self.table['MASS']
        else:
            x = self.table['X'][mask]
            y = self.table['Y'][mask]
            z = self.table['Z'][mask]
            w = self.table['MASS'][mask]
        coords = np.transpose(np.array([x, y, z]))
        edges: List
        masses, edges = np.histogramdd(coords, bins=bins, weights=w)
        d_bin = (edges[0][1] - edges[0][0]) / 2
        centers = np.array([i[:-1] for i in edges]) + d_bin
        self.d_bin = d_bin
#        print('d_bin = ', d_bin)
        masses = gaussian_filter(masses, 0.3)
        return masses, centers, d_bin

    def coordinates(self):
        coords = np.meshgrid(*self.centers)
        return np.asarray(coords)

    def distances(self, point):
        coords = self.coordinates()
        r = np.sqrt(np.sum([np.square(coords[i] - point[i]) for i in range(3)], 0))
        return r


    def potential(self, smooth: float = 0.0):
        g = cons.G
        u = np.zeros((self.n_bins,) * 3)
        for i, j in np.ndenumerate(u):
            s = (Ellipsis,) + i
            d = self.distances(self.coords[s])
            d[s] = self.d_bin / 10.0
            u[i] = -g * np.sum(self.masses / d)

        # if smooth != 0.0:
        #     u = gaussian_filter(u, smooth)

        index_min = np.where(u == u.min())
        coords_min = [float(self.coords[i][index_min]) for i in range(3)]

        print('Lowest potential coordinates: ({:.2f}, {:.2f}, {:.2f})'.format(*coords_min))
        print('Value for minimum potential = ', u.min())
        return u, np.asarray(coords_min)


if __name__ == '__main__':
    t_path = "/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/tables_arp142_v2"
    snapshot = 71
    n_bins = 12

    t = MyTable(table_path=t_path, snapshot=snapshot, galaxy=1, n_bins=n_bins)
    potential, coords_min = t.potential(smooth = 0.3)
    masses, centers, d_bin = t.binning(n_bins)

    positions1 = position_all_particles(t_path, snapshot=snapshot, galaxy=1)
    positions2 = position_all_particles(t_path, snapshot=snapshot, galaxy=2)
    gas_from_spiral = position(t_path, snapshot=snapshot, galaxy=2, p_type=0)

    plt.plot(positions1[0,:], positions1[1,:],',',label='Particles from eliptical galaxy')
    plt.plot(positions2[0,:], positions2[1,:],',',label='Particles from spiral galaxy')
    plt.plot(gas_from_spiral[0,:], gas_from_spiral[1,:],',',label='Gas from spiral galaxy')
    plt.title('snapshot ' + str(snapshot))
    plt.xlabel('X (KPc)')
    plt.ylabel('Y (Kpc)')

    plt.plot([coords_min[0] - d_bin, coords_min[0] + d_bin], [coords_min[1] - d_bin, coords_min[1] - d_bin], color='red')
    plt.plot([coords_min[0] - d_bin, coords_min[0] - d_bin], [coords_min[1] + d_bin, coords_min[1] - d_bin], color='red')
    plt.plot([coords_min[0] + d_bin, coords_min[0] + d_bin], [coords_min[1] - d_bin, coords_min[1] + d_bin], color='red')
    plt.plot([coords_min[0] + d_bin, coords_min[0] - d_bin], [coords_min[1] + d_bin, coords_min[1] + d_bin], color='red')

    plt.legend()
    plt.show()

    for i in range(10):
        mask = t.select_particles_cube(coords_min)
        t.masses, t.centers, d_bin = t.binning(n_bins, mask)
        t.coords = t.coordinates()
        potential, coords_min = t.potential(smooth = 0.3)
        plt.plot([coords_min[0] - d_bin, coords_min[0] + d_bin], [coords_min[1] - d_bin, coords_min[1] - d_bin], color='red')
        plt.plot([coords_min[0] - d_bin, coords_min[0] - d_bin], [coords_min[1] + d_bin, coords_min[1] - d_bin], color='red')
        plt.plot([coords_min[0] + d_bin, coords_min[0] + d_bin], [coords_min[1] - d_bin, coords_min[1] + d_bin], color='red')
        plt.plot([coords_min[0] + d_bin, coords_min[0] - d_bin], [coords_min[1] + d_bin, coords_min[1] + d_bin], color='red')
        plt.show()
        input()


    #fig = plt.figure()

    # for i in range(10):
    #     plt.clf()
    #     ax = fig.add_subplot(111)
    #     mask = t.select_particles_cube(coords_min)
    #     t.masses, t.centers = t.binning(10, mask)
    #     t.coords = t.coordinates()
    #     potential, coords_min = t.potential()
    #     print('coords = ', coords_min)
    #     ext = np.ravel([[t.coords[i].min(), t.coords[i].max()] for i in range(3)])
    #     print('ext = ', ext)
    #     ax.imshow(potential.sum(0), extent=ext[2:])
    #     plt.draw()
    #     input()


    # minimum_planes = np.where(potential == potential.min())
    #
    # fig, ax = plt.subplots(1, 3, sharey='row', figsize=(12, 4))
    # for j in range(3):
    #     idx = 3 * [slice(None)]
    #     idx[j] = int(minimum_planes[j])
    #     im = ax[j].imshow(potential[tuple(idx)], cmap='plasma_r', origin='lower')
    # plt.show()
