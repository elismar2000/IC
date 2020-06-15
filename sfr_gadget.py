from __future__ import print_function

import numpy as np
from pygadgetreader import *
from tables import *

class SFR:
    def __init__(self, table_mins, path_snapshot):
        self.table_mins = table_mins
        self.path = path_snapshot
        self.snapshot = str()
        self.previous_snapshot = str()
        self.snap_path = str()
        self.previous_snap_path = str()
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
        self.sm1 = np.array([])
        self.sm2 = np.array([])
        self.sm = np.array([])
        self.n1 = np.array([])
        self.n2 = np.array([])
        self.n = np.array([])


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
        Read positions of star paticles and their masses for
        each snapshot. To compute the masses of stars formed
        only during the last 50Myr, we have to mask all the star
        particles formed until the last snapshot
        '''

        if self.previous_snapshot == '0000':
            previous_ids = np.array([])

        else:
            previous_ids = readsnap(self.previous_snap_path, 'pid', 'star')

        current_ids = readsnap(self.snap_path, 'pid', 'star')

        new_born = np.isin(current_ids, previous_ids)

        mask = new_born == False

        masses = readsnap(self.snap_path, 'mass', 'star')
        pos = readsnap(self.snap_path, 'pos', 'star')

        self.masses = masses[mask]
        self.pos = pos[mask]


    def _distances(self, pos, min):
        '''
        Calculates the distance of the particles to the potential minima
        '''
        dist = np.sqrt(np.sum(np.square(pos[i] - min[i]) for i in range(3)))
        return dist[0]


    def sfr(self, radius, snapshot):
        '''
        Calculates the mass of gas inside a given radius from
        the potential minima of each galaxy for each snapshot
        '''
        self.snapshot = '%04d' % (snapshot,)
        self.previous_snapshot = '%04d' % (snapshot - 1,)

        self.snap_path = self.path + self.snapshot
        self.previous_snap_path = self.path + self.previous_snapshot

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

            self.sm1 = np.sum(self.masses[mask1])
            self.sm2 = np.sum(self.masses[mask2])

            self.n1 = sum(mask1)
            self.n2 = sum(mask2)

            print("mass of formed stars inside radius %g" %radius, "Kpc from coords_min1 = %g" %self.sm1)
            print("mass of formed stars inside radius %g" %radius, "Kpc from coords_min2 = %g" %self.sm2)

        else:
            coords_min = np.array([self.xmin, self.ymin, self.zmin])

            print("Computing particle distances from coords_min")
            dist = np.array([self._distances(self.pos[i], coords_min) for i in range(np.size(self.pos, 0))])

            mask = dist < radius

            self.sm = np.sum(self.masses[mask])

            self.n = sum(mask)

            print("mass of formed stars inside radius %g" %radius, "Kpc from coords_min = %g" %self.sm)


    def save_table_sfr(self, num_of_snaps):
        '''
        Write and save a hdf5 table with data from sfr
        The iteration must start from the second snapshot,
        i.e., snapshot_0001, due to the fact that the first snapshot
        have no stars yet formed.
        '''
        class Table(IsDescription):
            '''
            Class describing attributes for table in which to save
            sfr data
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

        h5file = open_file('galmer-like_sim_sfr.h5', mode='w', title='SFR table')
        group = h5file.create_group("/", 'sfr', 'Star formation rate')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        sfr = table.row

        for snapshot in range(1, num_of_snaps, 1):
            sfr['snapshot'] = snapshot
            print("This is snapshot ", snapshot)

            for radius in [1.0, 0.5, 0.1, 0.01]:
                self.sfr(radius, snapshot)
                print("Calculating sfr for radius ", radius, "Kpc")

                if self.two_mins:
                    if radius == 1.0:
                        sfr['gal1_r10'] = self.sm1
                        sfr['gal2_r10'] = self.sm2
                        sfr['gal_r10'] = np.nan

                    if radius == 0.5:
                        sfr['gal1_r05'] = self.sm1
                        sfr['gal2_r05'] = self.sm2
                        sfr['gal_r05'] = np.nan

                    if radius == 0.1:
                        sfr['gal1_r01'] = self.sm1
                        sfr['gal2_r01'] = self.sm2
                        sfr['gal_r01'] = np.nan

                    if radius == 0.01:
                        sfr['gal1_r001'] = self.sm1
                        sfr['gal2_r001'] = self.sm2
                        sfr['gal_r001'] = np.nan

                else:
                    if radius == 1.0:
                        sfr['gal1_r10'] = np.nan
                        sfr['gal2_r10'] = np.nan
                        sfr['gal_r10'] = self.sm

                    if radius == 0.5:
                        sfr['gal1_r05'] = np.nan
                        sfr['gal2_r05'] = np.nan
                        sfr['gal_r05'] = self.sm

                    if radius == 0.1:
                        sfr['gal1_r01'] = np.nan
                        sfr['gal2_r01'] = np.nan
                        sfr['gal_r01'] = self.sm

                    if radius == 0.01:
                        sfr['gal1_r001'] = np.nan
                        sfr['gal2_r001'] = np.nan
                        sfr['gal_r001'] = self.sm

            sfr.append()
        table.flush()


    def save_table_number(self, num_of_snaps):
        '''
        Write and save a hdf5 table with data from part_num
        '''
        class Table(IsDescription):
            '''
            Class describing attributes for table in which to save
            part_num data
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

        h5file = open_file('galmer-like_sim_part_num.h5', mode='w', title='Number of Particles Table')
        group = h5file.create_group("/", 'part_num', 'Number of Particles')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        part_num = table.row

        for snapshot in range(0, num_of_snaps, 1):
            part_num['snapshot'] = snapshot
            print("This is snapshot ", snapshot)

            for radius in [1.0, 0.5, 0.1, 0.01]:
                self.inflow(radius, snapshot)
                print("Calculating inflow for radius ", radius, "Kpc")

                if self.two_mins:
                    if radius == 1.0:
                        part_num['gal1_r10'] = self.n1
                        part_num['gal2_r10'] = self.n2
                        part_num['gal_r10'] = np.nan

                    if radius == 0.5:
                        part_num['gal1_r05'] = self.n1
                        part_num['gal2_r05'] = self.n2
                        part_num['gal_r05'] = np.nan

                    if radius == 0.1:
                        part_num['gal1_r01'] = self.n1
                        part_num['gal2_r01'] = self.n2
                        part_num['gal_r01'] = np.nan

                    if radius == 0.01:
                        part_num['gal1_r001'] = self.n1
                        part_num['gal2_r001'] = self.n2
                        part_num['gal_r001'] = np.nan

                else:
                    if radius == 1.0:
                        part_num['gal1_r10'] = np.nan
                        part_num['gal2_r10'] = np.nan
                        part_num['gal_r10'] = self.n

                    if radius == 0.5:
                        part_num['gal1_r05'] = np.nan
                        part_num['gal2_r05'] = np.nan
                        part_num['gal_r05'] = self.n

                    if radius == 0.1:
                        part_num['gal1_r01'] = np.nan
                        part_num['gal2_r01'] = np.nan
                        part_num['gal_r01'] = self.n

                    if radius == 0.01:
                        part_num['gal1_r001'] = np.nan
                        part_num['gal2_r001'] = np.nan
                        part_num['gal_r001'] = self.n

            part_num.append()
        table.flush()


if __name__ == '__main__':
    table_mins = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins_Gadget/galmer-like_sim_test2_minima.h5'
    path_snapshot = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/snapshot_'

    s = SFR(table_mins, path_snapshot)
    s.save_table_sfr(num_of_snaps=71)
