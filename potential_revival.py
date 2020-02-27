import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import scipy.constants as cons


class MyPotential:
    def __init__(self, table_path: str, bins: int = 25):
        self.table_path = table_path
        self.bins = bins
        self.table = Table.read(self.table_path)
        self.x = self.table['X']
        self.y = self.table['Y']
        self.m = self.table['MASS']
        self.mask = None
        self.masses, self.edges = self.binning()
        self.center, self.d_binx, self.d_biny = self.centers()
        self.coords = self.coordinates()


    def select_particles_cube(self, coords):
        widthx = self.d_binx * 1.5
        widthy = self.d_biny * 1.5
        mask = (self.x >= (coords[0] - widthx)) & (self.x <= (coords[0] + widthx))\
        & (self.y >= (coords[1] - widthy)) & (self.y <= (coords[1] + widthy))
        self.mask = mask


    def binning(self):
        if self.mask is None:
            x = self.x
            y = self.y
            m = self.m
        else:
            x = self.x[self.mask]
            y = self.y[self.mask]
            m = self.m[self.mask]

        coords = np.transpose(np.array([x, y]))
        masses, edges = np.histogramdd(coords, bins=self.bins, weights=m)
        edges = np.asarray(edges)
        return masses, edges


    def centers(self):
        x = self.x
        y = self.y
        self.masses, self.edges = self.binning()
        d_binx = (x.max() - x.min()) / self.bins
        d_biny = (y.max() - y.min()) / self.bins
        a = self.edges[0] + d_binx / 2
        b = self.edges[1] + d_biny / 2
        center_x = a[:-1]
        center_y = b[:-1]
        centers = np.array([center_x, center_y])

        # print('x = ', x)
        # print('y = ', y)
        # print('edges = ', edges)
        return centers, d_binx, d_biny


    def coordinates(self):
        self.center, self.d_binx, self.d_biny = self.centers()
        coords = np.meshgrid(*self.center)
        return np.asarray(coords)


    def distances(self, point):
        self.coords = self.coordinates()
        r = np.sqrt(np.sum([np.square(self.coords[i] - point[i]) for i in range(2)], 0))
        return r


    def potential(self):
        g = cons.G
        u = np.zeros((self.bins,) * 2)
        dx = self.d_binx / 2
        dy = self.d_biny / 2

        for i, j in np.ndenumerate(u):
            s = (Ellipsis,) + i
            d = self.distances(self.coords[s])
            d[s] = np.sqrt(dx**2 + dy**2) / 2.0
            u[i] = -g * np.sum(self.masses / d)

        index_min = np.where(u == u.min())
        coords_min = [float(self.coords[i][index_min]) for i in range(2)]
        import pdb; pdb.set_trace()
        print('Lowest potential coordinates: ({:.2f}, {:.2f})'.format(*coords_min))
        print('Value for minimum potential = ', u.min())

        return u, coords_min


if __name__ == '__main__':

    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/gas_velocity/fake_galaxy.fits'
    bins = 10
    t = MyPotential(table_path=t_path, bins=bins)
    table = Table.read(t_path)

    potential, coords_min = t.potential()
    dx = t.d_binx / 2
    dy = t.d_biny / 2

    ext = [t.edges[0, 0], t.edges[0, 10], t.edges[1, 0], t.edges[1, 10]]
    plt.imshow(potential, cmap='plasma_r', origin='lower', extent=ext)
    plt.colorbar()
    plt.show()

    # plt.hist2d(table['X'], table['Y'], bins=100, weights=table['MASS'])
    # plt.plot(coords_min[0], coords_min[1], 'X')
    # plt.plot([coords_min[0] - dx, coords_min[0] + dx], [coords_min[1] - dy, coords_min[1] - dy], color='red')
    # plt.plot([coords_min[0] - dx, coords_min[0] - dx], [coords_min[1] + dy, coords_min[1] - dy], color='red')
    # plt.plot([coords_min[0] + dx, coords_min[0] + dx], [coords_min[1] - dy, coords_min[1] + dy], color='red')
    # plt.plot([coords_min[0] + dx, coords_min[0] - dx], [coords_min[1] + dy, coords_min[1] + dy], color='red')
    # plt.show()
    #
    #
    for i in range(10):
        t.select_particles_cube(coords_min)
        potential, coords_min = t.potential()
        dx = t.d_binx / 2
        dy = t.d_biny / 2
        plt.plot(coords_min[0], coords_min[1], 'X')
        plt.plot([coords_min[0] - dx, coords_min[0] + dx], [coords_min[1] - dy, coords_min[1] - dy], color='red')
        plt.plot([coords_min[0] - dx, coords_min[0] - dx], [coords_min[1] + dy, coords_min[1] - dy], color='red')
        plt.plot([coords_min[0] + dx, coords_min[0] + dx], [coords_min[1] - dy, coords_min[1] + dy], color='red')
        plt.plot([coords_min[0] + dx, coords_min[0] - dx], [coords_min[1] + dy, coords_min[1] + dy], color='red')
        plt.draw()
        input()

    # fig = plt.figure()
    #
    # for i in range(10):
    #     plt.clf()
    #     ax = fig.add_subplot(111)
    #     mask = t.select_particles_cube(coords_min)
    #     t.masses, t.edges = t.binning(mask)
    #     t.centers, t.d_binx, t.d_biny = t.centers()
    #     t.coords = t.coordinates()
    #     potential, coords_min = t.potential()
    #     ext = np.ravel([[t.coords[i].min(), t.coords[i].max()] for i in range(2)])
    #     print('ext = ', ext)
    #     ax.imshow(potential.sum(0), extent=ext[2:])
    #     plt.draw()
    #     input()
