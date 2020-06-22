import astropy.table
import numpy as np
import pickle
from tables import *

import os

class SFR:
    def __init__(self, t_path: str, table_mins: str):
        self.t_path = t_path
        self.table_mins = table_mins
        self.snapshot = 0
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.sm1 = np.array([])
        self.sm2 = np.array([])
        self.sm = np.array([])
        self.previous_sm1 = np.array([])
        self.previous_sm2 = np.array([])
        self.previous_sm = np.array([])


    def read_table(self):
        '''
        Read mass of recently formed stars
        and their x, y, z positions for each snapshot
        '''
        t = astropy.table.Table.read(self.t_path, self.snapshot)

        gas = t['P_TYPE'] == 0
        new_born = t['M_GAS'][gas] < (0.05 * t['MASS'][gas])
        #mudar aqui

        self.x = t['X'][gas][new_born]
        self.y = t['Y'][gas][new_born]
        self.z = t['Z'][gas][new_born]

        self.mass = 0.95 * t['MASS'][gas][new_born]


    def sfr(self, radius: float = 1.0, merge: int = 17):
        '''
        Computes mass of formed stars inside a given radius
        from the potential minima.

        Parameters
        -----------
        radius : float
            Radius from the potential minima within which
            we sum up the masses of star particles

        merge : int
            The last snapshot with galaxies still not coalesced
        '''
        self.read_table()

        r = np.array([self.x, self.y, self.z])
        print('In snapshot ', self.snapshot, 'r = ', r)

        with open(self.table_mins, 'rb') as f:
            table = pickle.load(f)

        if self.snapshot <= merge:
            center1 = table[self.snapshot - 1, 1]
            center2 = table[self.snapshot - 1, 2]

            dist1 = np.array([self._distances(r[:, i], center1) for i in range(len(r[0]))])
            dist2 = np.array([self._distances(r[:, i], center2) for i in range(len(r[0]))])

            mask1 = dist1 < radius
            mask2 = dist2 < radius

            if self.snapshot == 2:
                self.sm1 = np.sum(self.mass[mask1])
                self.sm2 = np.sum(self.mass[mask2])

                self.previous_sm1 = self.sm1
                self.previous_sm2 = self.sm2

            else:
                self.sm1 = np.sum(self.mass[mask1]) - self.previous_sm1
                self.sm2 = np.sum(self.mass[mask2]) - self.previous_sm2

                self.previous_sm1 = self.sm1
                self.previous_sm2 = self.sm2

        if self.snapshot > merge:
            center = table[self.snapshot - 1, 3]

            dist = np.array([self._distances(r[:, i], center) for i in range(len(r[0]))])

            mask = dist < radius

            if self.snapshot == merge + 1:
                self.sm = np.sum(self.mass[mask]) - (self.previous_sm1 + self.previous_sm2)

                self.previous_sm = self.sm

            else:
                self.sm = np.sum(self.mass[mask]) - self.previous_sm


    def _distances(self, pos, min):
        '''
        Calculates the distance of the particles to the potential minima
        '''
        dist = np.sqrt(np.sum(np.square(pos[i] - min[i]) for i in range(3)))
        return dist


    def save_table(self, num_of_snaps):
        '''
        Save a hdf5 table with the data from sfr
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

        h5file = open_file('orbit1_galmer_sim_sfr.h5', mode='w', title='SFR table')
        group = h5file.create_group("/", 'sfr', 'Star formation rate')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        sfr = table.row

        for snapshot in range(2, num_of_snaps + 1, 1):
            self.snapshot = snapshot

            sfr['snapshot'] = snapshot
            print("This is snapshot ", snapshot)

            for radius in [1.0, 0.5, 0.1, 0.01]:
                merge = 17
                self.sfr(radius, merge)
                print("Calculating sfr for radius ", radius, "Kpc")

                if self.snapshot <= merge:
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


if __name__ == '__main__':
    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit1'
    orbit1_mins = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/mins_GalMer/orbit1_mins.txt'

    s = SFR(t_path, orbit1_mins)
    s.save_table(num_of_snaps=71)
