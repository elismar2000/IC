import matplotlib.pyplot as plt
import numpy as np
from pygadgetreader import *
from tables import *

class Potential:
    def __init__(self, path):
        self.path = path
        self.snapshot = str()
        self.snap_path = str()
        self.pos = np.array([])
        self.pot = np.array([])
        self.pos_stars = np.array([])
        self.pot_stars = np.array([])
        self.pid_stars = np.array([])
        self.coords_min1 = np.array([])
        self.coords_min2 = np.array([])
        self.coords_min = np.array([])
        self.two_min = True
        self.id_min1 = np.array([])
        self.id_min2 = np.array([])


    def positions(self):
        '''
        Read particles' positions
        '''
        pos_halo = readsnap(self.snap_path, 'pos', 'dm')
        pos_disk = readsnap(self.snap_path, 'pos', 'disk')
        pos_gas = readsnap(self.snap_path, 'pos', 'gas')
        pos_bulge = readsnap(self.snap_path, 'pos', 'bulge')

        if self.stars():
            pos_stars = readsnap(self.snap_path, 'pos', 'star')
            self.pos = np.concatenate((pos_halo, pos_disk, pos_gas, pos_bulge, pos_stars))

        else:
            self.pos = np.concatenate((pos_halo, pos_disk, pos_gas, pos_bulge))


    def ids(self):
        '''
        Read particles' ids
        '''
        id_halo = readsnap(self.snap_path, 'pid', 'dm')
        id_disk = readsnap(self.snap_path, 'pid', 'disk')
        id_gas = readsnap(self.snap_path, 'pid', 'gas')
        id_bulge = readsnap(self.snap_path, 'pid', 'bulge')

        if self.stars():
            id_stars = readsnap(self.snap_path, 'pid', 'star')
            self.pid = np.concatenate((id_halo, id_disk, id_gas, id_bulge, id_stars))

        else:
            self.pid = np.concatenate((id_halo, id_disk, id_gas, id_bulge))


    def potential(self, beta, gama):
        '''
        Find gravitational potential minima for each galaxy in each snapshot
        of the simulation. First, just determines what is the particle with
        the lower value of gravitational potential (galaxy 1) and then masks a region of
        radius 'beta' around it, and finds the particle with lower gravitational
        potential among the rest of the particles (galaxy 2). In order to ensure
        that each potential minima in each snapshot will be in the right column of
        the table created at the end of the code, in each iteration the ids of the
        particles in a region inside a radius 'beta' are recorded to be used during
        the next iteration. If the ids of the majority of the particles inside that
        same radius in the forthcoming iteration are equal to the ids recorded in the
        previous iteration, them they are representing the same galaxy. If not, they
        are representing the other galaxy.

        Parameters
        ----------
        beta : float
            The radius used to create a mask around each potential minima, applied
            both to find the potential minima of the other galaxy and to record the ids
            of the particles inside a sphere with such a radius.
            beta = 1.0 (Kpc) usually works

        gama : float
            if pot2 * gama < pot1, then it's considered to be the potential minima
            of the second galaxy. If not, then there is just one minimum in the system.
            pot1 and pot2 are the values of the gravitational potential minima of each galaxy.
            gama = 0.75 usually works
        '''

        pot_halo = readsnap(self.snap_path, 'pot', 'dm')
        pot_disk = readsnap(self.snap_path, 'pot', 'disk')
        pot_gas = readsnap(self.snap_path, 'pot', 'gas')
        pot_bulge = readsnap(self.snap_path, 'pot', 'bulge')

        if self.stars():
            pot_stars = readsnap(self.snap_path, 'pot', 'star')
            self.pot = np.concatenate((pot_halo, pot_disk, pot_gas, pot_bulge, pot_stars))

        else:
            self.pot = np.concatenate((pot_halo, pot_disk, pot_gas, pot_bulge))

        self.positions()
        self.ids()

        if self.snapshot == 'snapshot_0000':
            pot1 = self.pot.min()
            print('pot1: ', pot1)
            index_min1 = np.where(self.pot == pot1)[0][0]
            self.coords_min1 = self.pos[index_min1]

            dist1 = np.array([self._distances(self.pos[i], self.coords_min1) for i in range(len(self.pos))])
            mask_out = dist1 > beta
            mask_in1 = dist1 < beta

            pot2 = self.pot[mask_out].min()
            print('pot2: ', pot2)
            index_min2 = np.where(self.pot == pot2)[0][0]
            self.coords_min2 = self.pos[index_min2]
            print('There are two min')

            dist2 = np.array([self._distances(self.pos[i], self.coords_min2) for i in range(len(self.pos))])
            mask_in2 = dist2 < beta

            self.id_min1 = self.pid[mask_in1]
            self.id_min2 = self.pid[mask_in2]

        else:
            pot1 = self.pot.min()
            print('pot1: ', pot1)
            index_min1 = np.where(self.pot == pot1)[0][0]
            coords_min1 = self.pos[index_min1]

            dist1 = np.array([self._distances(self.pos[i], coords_min1) for i in range(len(self.pos))])
            mask_out1 = dist1 > beta
            mask_in1 = dist1 < beta

            pot2 = self.pot[mask_out1].min()
            print('pot2: ', pot2)

            self.two_min = pot2 < pot1 * gama
            if self.two_min:
                index_min2 = np.where(self.pot == pot2)[0][0]
                coords_min2 = self.pos[index_min2]
                print('There are two min')

                dist2 = np.array([self._distances(self.pos[i], coords_min2) for i in range(len(self.pos))])
                mask_in2 = dist2 < beta

                current_id_min1 = self.pid[mask_in1]
                current_id_min2 = self.pid[mask_in2]

                is1 = np.isin(current_id_min1, self.id_min1)
                is2 = np.isin(current_id_min1, self.id_min2)

                print('sum is1, is2 = ', [sum(is1), sum(is2)])
                if sum(is1) > sum(is2):
                    self.coords_min1 = coords_min1
                    self.coords_min2 = coords_min2

                    self.id_min1 = current_id_min1
                    self.id_min2 = current_id_min2
                    print('1 = 1, 2 = 2')

                else:
                    self.coords_min1 = coords_min2
                    self.coords_min2 = coords_min1

                    self.id_min1 = current_id_min2
                    self.id_min2 = current_id_min1
                    print('1 = 2, 2 = 1')

            else:
                self.coords_min = coords_min1
                print('There is only one min')

        #self.plot()


    def plot(self):
        '''
        Plot positions of potential minima
        '''
        plt.figure()
        plt.plot(self.pos[:, 0], self.pos[:, 1], ',')

        if self.two_min:
            plt.plot(self.coords_min1[0], self.coords_min1[1], 'y*')
            plt.plot(self.coords_min2[0], self.coords_min2[1], 'rX')

        else:
            plt.plot(self.coords_min[0], self.coords_min[1], 'rX')

        plt.show()


    def _distances(self, position, point):
        '''
        Evaluates distances of particles from a given point
        '''
        distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
        return distances


    def stars(self):
        '''
        Read snapshot atributes for formed stars
        '''
        header = readheader(self.snap_path, 'header')
        nstars = header['nstar']
        print('nstars = ', nstars)

        if nstars != 0:
            return True
        else:
            return False


    def run(self, num_of_snaps, interval):
        '''
        perform the iteration over all snapshots
        '''
        class Table(IsDescription):
            '''
            Class describing attributes for table in which to save
            potential min coordinates
            '''
            snapshot = Int32Col()
            xmin1 = Float32Col()
            ymin1 = Float32Col()
            zmin1 = Float32Col()
            xmin2 = Float32Col()
            ymin2 = Float32Col()
            zmin2 = Float32Col()
            xmin = Float32Col()
            ymin = Float32Col()
            zmin = Float32Col()

        h5file = open_file('galmer-like_sim_minima2.h5', mode='w', title='Min table')
        group = h5file.create_group("/", 'potential', 'Gravitational potential minima')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        pot = table.row

        for snapshot in range(0, num_of_snaps, interval):
            self.snapshot = 'snapshot_' + '%04d' % (snapshot,)
            self.snap_path = self.path + self.snapshot
            print('snapshot: ', self.snapshot)
            self.potential(beta=1.0, gama=0.75)

            if p.two_min:
                pot['snapshot'] = snapshot
                pot['xmin1'] = self.coords_min1[0]
                pot['ymin1'] = self.coords_min1[1]
                pot['zmin1'] = self.coords_min1[2]
                pot['xmin2'] = self.coords_min2[0]
                pot['ymin2'] = self.coords_min2[1]
                pot['zmin2'] = self.coords_min2[2]
                pot['xmin'] = np.nan
                pot['ymin'] = np.nan
                pot['zmin'] = np.nan
                pot.append()

            else:
                pot['snapshot'] = snapshot
                pot['xmin1'] = np.nan
                pot['ymin1'] = np.nan
                pot['zmin1'] = np.nan
                pot['xmin2'] = np.nan
                pot['ymin2'] = np.nan
                pot['zmin2'] = np.nan
                pot['xmin'] = self.coords_min[0]
                pot['ymin'] = self.coords_min[1]
                pot['zmin'] = self.coords_min[2]
                pot.append()

        table.flush()


if __name__ == '__main__':
    path = '/home/elismar/Documentos/Fisica/IC/Gadget3/simulation_galmer-like_test2/'

    p = Potential(path)
    p.run(71, 1)
