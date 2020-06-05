from __future__ import print_function

import numpy as np
from pygadgetreader import *
from tables import *

class Inflow:
    def __init__(self, table_mins, path_snapshot):
        self.table_mins = table_mins
        self.path = path_snapshot
        self.snapshot = str()
        self.snap_path = str()
        self.xmin1 = np.array([])
        self.ymin1 = np.array([])
        self.zmin1 = np.array([])
        self.xmin2 = np.array([])
        self.ymin2 = np.array([])
        self.zmin2 = np.array([])
        self.xmin = np.array([])
        self.ymin = np.array([])
        self.zmin = np.array([])
        self.two_mins = True
        self.pos = np.array([])
        self.masses = np.array([])
        self.gm1 = np.array([])
        self.gm2 = np.array([])
        self.gm = np.array([])


    def coords_min(self, snap):
        '''
        Read table of gravitational potential minima
        and return their positions for a given snapshot
        '''
        h5file = open_file(self.table_mins, mode='r')
        table = h5file.root.potential.readout

        self.xmin1 = np.array([j['xmin1'] for j in table.where('snapshot == snap')])
        self.ymin1 = np.array([j['ymin1'] for j in table.where('snapshot == snap')])
        self.zmin1 = np.array([j['zmin1'] for j in table.where('snapshot == snap')])

        self.xmin2 = np.array([j['xmin2'] for j in table.where('snapshot == snap')])
        self.ymin2 = np.array([j['ymin2'] for j in table.where('snapshot == snap')])
        self.zmin2 = np.array([j['zmin2'] for j in table.where('snapshot == snap')])

        self.xmin = np.array([j['xmin'] for j in table.where('snapshot == snap')])
        self.ymin = np.array([j['ymin'] for j in table.where('snapshot == snap')])
        self.zmin = np.array([j['zmin'] for j in table.where('snapshot == snap')])

        if np.isnan(self.xmin):
            self.two_mins = True
            print("There are two galaxies in this snapshot")

        if np.isnan(self.xmin1):
            self.two_mins = False
            print("There is just one galaxy is this snapshot")

        h5file.close()


    def read_snapshots(self):
        '''
        Read positions of gas paticles and their masses for
        each snapshot
        '''
        self.pos = readsnap(self.snap_path, 'pos', 'gas')
        self.masses = readsnap(self.snap_path, 'mass', 'gas')


    def _distances(self, pos, min):
        '''
        Calculates the distance of the particles to the potential minima
        '''
        dist = np.sqrt(np.sum(np.square(pos[i] - min[i]) for i in range(3)))
        return dist[0]


    def inflow(self, radius, snapshot):
        '''
        Calculates the mass of gas inside a given radius from
        the potential minima of each galaxy for each snapshot
        '''
        self.snapshot = '%04d' % (snapshot,)
        self.snap_path = self.path + self.snapshot

        self.read_snapshots()
        self.coords_min(snapshot)

        if self.two_mins:
            coords_min1 = np.array([self.xmin1, self.ymin1, self.zmin1])
            coords_min2 = np.array([self.xmin2, self.ymin2, self.zmin2])

            print("Computing particle distances from coords_min1")
            dist1 = np.array([self._distances(self.pos[i], coords_min1) for i in range(np.size(self.pos, 0))])

            print("Computing particle distances from coords_min2")
            dist2 = np.array([self._distances(self.pos[i], coords_min2) for i in range(np.size(self.pos, 0))])

            mask1 = dist1 < radius
            mask2 = dist2 < radius

            self.gm1 = np.sum(self.masses[mask1])
            self.gm2 = np.sum(self.masses[mask2])

            print("gas mass inside radius %g" %radius, "Kpc from coords_min1 = %g" %self.gm1)
            print("gas mass inside radius %g" %radius, "Kpc from coords_min2 = %g" %self.gm2)

        else:
            coords_min = np.array([self.xmin, self.ymin, self.zmin])

            print("Computing particle distances from coords_min")
            dist = np.array([self._distances(self.pos[i], coords_min) for i in range(np.size(self.pos, 0))])

            mask = dist < radius

            self.gm = np.sum(self.masses[mask])

            print("gas mass inside radius %g" %radius, "Kpc from coords_min = %g" %self.gm)


    def save_table(self, num_of_snaps):
        '''
        Write and save a hdf5 table with data from inflow
        '''
        class Table(IsDescription):
            '''
            Class describing attributes for table in which to save
            inflow data
            '''
            snapshot = Int32Col()
            gal1_r10 = Float32Col()
            gal1_r05 = Float32Col()
            gal1_r01 = Float32Col()
            gal1_r001 = Float32Col()
            gal2_r10 = Float32Col()
            gal2_r05 = Float32Col()
            gal2_r01 = Float32Col()
            gal2_r001 = Float32Col()
            gal_r10 = Float32Col()
            gal_r05 = Float32Col()
            gal_r01 = Float32Col()
            gal_r001 = Float32Col()

        h5file = open_file('galmer-like_sim_inflow.h5', mode='w', title='Inflow table')
        group = h5file.create_group("/", 'inflow', 'Inflow of Gas')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        inflow = table.row

        for snapshot in range(0, num_of_snaps, 1):
            inflow['snapshot'] = snapshot
            print("This is snapshot ", snapshot)

            for radius in [1.0, 0.5, 0.1, 0.01]:
                self.inflow(radius, snapshot)
                print("Calculating inflow for radius ", radius, "Kpc")

                if self.two_mins:
                    if radius == 1.0:
                        inflow['gal1_r10'] = self.gm1
                        inflow['gal2_r10'] = self.gm2
                        inflow['gal_r10'] = np.nan

                    if radius == 0.5:
                        inflow['gal1_r05'] = self.gm1
                        inflow['gal2_r05'] = self.gm2
                        inflow['gal_r05'] = np.nan

                    if radius == 0.1:
                        inflow['gal1_r01'] = self.gm1
                        inflow['gal2_r01'] = self.gm2
                        inflow['gal_r01'] = np.nan

                    if radius == 0.01:
                        inflow['gal1_r001'] = self.gm1
                        inflow['gal2_r001'] = self.gm2
                        inflow['gal_r001'] = np.nan

                else:
                    if radius == 1.0:
                        inflow['gal1_r10'] = np.nan
                        inflow['gal2_r10'] = np.nan
                        inflow['gal_r10'] = self.gm

                    if radius == 0.5:
                        inflow['gal1_r05'] = np.nan
                        inflow['gal2_r05'] = np.nan
                        inflow['gal_r05'] = self.gm

                    if radius == 0.1:
                        inflow['gal1_r01'] = np.nan
                        inflow['gal2_r01'] = np.nan
                        inflow['gal_r01'] = self.gm

                    if radius == 0.01:
                        inflow['gal1_r001'] = np.nan
                        inflow['gal2_r001'] = np.nan
                        inflow['gal_r001'] = self.gm

            inflow.append()
        table.flush()


if __name__ == '__main__':
    table_mins = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins_Gadget/galmer-like_sim_test2_minima.h5'
    path_snapshot = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_'

    i = Inflow(table_mins, path_snapshot)
    i.save_table(num_of_snaps=71)
