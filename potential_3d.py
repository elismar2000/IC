import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.table import Table
from relative_minima import find_minima
from scipy.optimize import minimize
from scipy.interpolate import interp1d

class Potential:
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.mass = np.array([])
        self.mask = np.array([])
        self.dx = 0.0
        self.dy = 0.0
        self.dz = 0.0

    def read_table(self, table_path, snapshot: str = None):
        if snapshot != None:
            t = Table.read(table_path, snapshot)
            self.x = t['X']
            self.y = t['Y']
            self.z = t['Z']
            self.mass = t['MASS']

        if snapshot == None:
            t = Table.read(table_path)
            self.x = t['X']
            self.y = t['Y']
            self.z = t['Z']
            self.mass = t['MASS']


    def evaluate_potential(self, n_bins: int = 10, smooth: float = 1.3):
        masses, edges = np.histogramdd(np.column_stack([self.x.data, self.y.data, self.z.data]), bins=n_bins, weights=self.mass)
        n, edges2 = np.histogramdd(np.column_stack([self.x.data, self.y.data, self.z.data]), bins=n_bins)
        edges = np.asarray(edges)

        dx = (edges[0][-1] - edges[0][0]) / n_bins
        dy = (edges[1][-1] - edges[1][0]) / n_bins
        dz = (edges[2][-1] - edges[2][0]) / n_bins
        self.dx = dx
        self.dy = dy
        self.dz = dz
        center = np.array([edges[0][:-1] + dx / 2, edges[1][:-1] + dy / 2, edges[2][:-1] + dz / 2])
        coords = np.meshgrid(*center)

        potential = np.zeros_like(masses)
        for i, u in np.ndenumerate(potential):
            distance = np.sqrt(np.sum([np.square(coords[j] - coords[j][i]) for j in range(3)], 0))
            distance[i] = dx / 4.0
            potential += -(masses / distance)

        if smooth != 0.0:
            for i in range(np.size(potential, axis=2)):
                potential[:, :, i] = convolve(potential[:, :, i], Gaussian2DKernel(smooth))

        return potential, center, n


    def select(self, coords, width: float = 2.0):
        wx = self.dx * width
        wy = self.dy * width
        wz = self.dz * width

        mask = (np.abs(self.x - coords[0]) <= wx) & (np.abs(self.y - coords[1]) <= wy) & (np.abs(self.z - coords[2]) <= wz)
        self.x = self.x[mask]
        self.y = self.y[mask]
        self.z = self.z[mask]
        self.mass = self.mass[mask]


    def pot3d(self, t_path, snapshot: int = None, n_bins: int = 20, smooth: float = 1.3, width: float = 1.5):
        self.read_table(t_path, snapshot)
        n_min = 200.0
        while n_min > 100.0:
            p, c, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)

            index_min = np.asarray(find_minima(p))
            mins = p[index_min[:, 0], index_min[:, 1], index_min[:, 2]]
            print('mins: ', mins)
            i = np.where(mins < mins.min() / 1.5)

            j = []
            for _ in range(np.size(i)):
                j.append(np.where(p == mins[i[0][_]]))

            if np.size(j, axis=0) == 1:
                coords_min = np.array([float(c[_][j[0][_]]) for _ in range(3)])
                self.select(coords_min, width=width)
                n_min = float(n[j[0]])
                print('Still finding minima')

            if np.size(j, axis=0) == 2:
                print('There was two minima')
                dx, dy, dz = self.dx, self.dy, self.dz
                coords_min1 = np.array([float(c[_][j[0][_]]) for _ in range(3)])
                coords_min1 = self.converge(coords_min1, t_path, snapshot, n_bins=n_bins, smooth=smooth, width=width)

                self.dx, self.dy, self.dz = dx, dy, dz
                coords_min2 = np.array([float(c[_][j[1][_]]) for _ in range(3)])
                coords_min2 = self.converge(coords_min2, t_path, snapshot, n_bins=n_bins, smooth=smooth, width=width)

                coords_min = np.array([coords_min1, coords_min2])

                break

            if np.size(j, axis=0) > 2:
                print('There is one minimum')
                i = np.where(p == p.min())
                coords_min = np.array([float(c[_][i[_]]) for _ in range(3)])
                self.select(coords_min, width=width)
                n_min = float(n[i])
                if n_min <= 100.0:
                    coords_min = np.array([self.spline(c, p, coords_min, i, dim='x'),
                                           self.spline(c, p, coords_min, i, dim='y'),
                                           self.spline(c, p, coords_min, i, dim='z')])
                    break

        return coords_min


    def converge(self, coords_min, t_path, snapshot, n_bins: int = 10, smooth: float = 1.3, width: float = 1.5):
        self.read_table(t_path, snapshot)
        self.select(coords_min, width=width)
        n_min = 200.0
        while n_min > 100.0:
            p, c, n = self.evaluate_potential(n_bins=n_bins, smooth=smooth)
            index_min = np.where(p == p.min())
            coords_min = np.array([float(c[i][index_min[i]]) for i in range(3)])
            self.select(coords_min, width=width)

            n_min = n[index_min]

        coords_min = np.array([self.spline(c, p, coords_min, index_min, dim='x'),
                               self.spline(c, p, coords_min, index_min, dim='y'),
                               self.spline(c, p, coords_min, index_min, dim='z')])

        return coords_min


    def spline(self, centers, pot, coords_min, index, dim, plot: bool = False):
        if dim == 'x':
            pot = pot[:, int(index[1]), int(index[2])]
            c = centers[0]
            x0 = coords_min[0]

            f = interp1d(c, pot, kind='cubic', fill_value='extrapolate')

            m = minimize(f, [x0], bounds=[[c.min(), c.max()]])
            min = float(m.x)

            if plot == True:
                plt.figure()
                cnew = np.linspace(c.min(), c.max(), num=100)
                plt.plot(c, pot, 'k--', label='Potencial discreto (x)')
                plt.plot(cnew, f(cnew), 'r-', label='Spline cubico')
                plt.plot(m.x, m.fun, 'X', label='Ponto de mínimo')
                plt.title('Ajuste do potencial em x')
                plt.legend()
                plt.show()

        if dim == 'y':
            pot = pot[int(index[0]), :, int(index[2])]
            c = centers[1]
            x0 = coords_min[1]

            f = interp1d(c, pot, kind='cubic', fill_value='extrapolate')

            m = minimize(f, [x0], bounds=[[c.min(), c.max()]])
            min = float(m.x)

            if plot == True:
                plt.figure()
                cnew = np.linspace(c.min(), c.max(), num=100)
                plt.plot(c, pot, 'k--', label='Potencial discreto (y)')
                plt.plot(cnew, f(cnew), 'r-', label='Spline cubico')
                plt.plot(m.x, m.fun, 'X', label='Ponto de mínimo')
                plt.title('Ajuste do potencial em y')
                plt.legend()
                plt.show()

        if dim == 'z':
            pot = pot[int(index[0]), int(index[1]), :]
            c = centers[2]
            x0 = coords_min[2]

            f = interp1d(c, pot, kind='cubic', fill_value='extrapolate')

            m = minimize(f, [x0], bounds=[[c.min(), c.max()]])
            min = float(m.x)

            if plot == True:
                plt.figure()
                cnew = np.linspace(c.min(), c.max(), num=100)
                plt.plot(c, pot, 'k--', label='Potencial discreto (z)')
                plt.plot(cnew, f(cnew), 'r-', label='Spline cubico')
                plt.plot(m.x, m.fun, 'X', label='Ponto de mínimo')
                plt.title('Ajuste do potencial em z')
                plt.legend()
                plt.show()

        return min

if __name__ == '__main__':

    p = Potential()

    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit1'
    snapshot = 1
    smooth = 1.3
    width = 1.0
    n_bins = 10
    coords_min = p.pot3d(t_path, snapshot, n_bins=n_bins, smooth=smooth, width=width)

    loop = False

    if loop == True:
        if coords_min.ndim == 1:
            table = np.array([snapshot, np.nan, np.nan, coords_min])
        if coords_min.ndim == 2:
            table = np.array([snapshot, coords_min[0], coords_min[1], np.nan])
        for j in range(snapshot + 1, 72):
            print('snapshot = ', j)
            coords_min = p.pot3d(t_path, snapshot=j, n_bins=n_bins, smooth=smooth, width=width)

            if coords_min.ndim == 2:
                if table.ndim == 1:
                    d1 = np.sqrt(np.sum([np.square(coords_min[0, i] - table[1][i]) for i in range(3)], 0))
                    d2 = np.sqrt(np.sum([np.square(coords_min[0, i] - table[2][i]) for i in range(3)], 0))
                if table.ndim >= 2:
                    d1 = np.sqrt(np.sum([np.square(coords_min[0, i] - table[j - (snapshot + 1), 1][i]) for i in range(3)], 0))
                    d2 = np.sqrt(np.sum([np.square(coords_min[0, i] - table[j - (snapshot + 1), 2][i]) for i in range(3)], 0))
                if d1 < d2:
                    table = np.vstack((table, np.array([j, coords_min[0], coords_min[1], np.nan])))
                if d1 > d2:
                    table = np.vstack((table, np.array([j, coords_min[1], coords_min[0], np.nan])))

            if coords_min.ndim == 1:
                table = np.vstack((table, np.array([j, np.nan, np.nan, coords_min])))

#########################################################################################################################

    if loop == False:
        from inflow import Inflow
        i = Inflow(t_path, snapshot)
        i.read_table()

        fig, axs = plt.subplots(1, 3, figsize=(20, 20))

        axs[0].plot(i.x, i.y, ',')
        if coords_min.shape == (3,):
            axs[0].plot(coords_min[0], coords_min[1], 'ko')
        if coords_min.shape == (2, 3):
            axs[0].plot(coords_min[0, 0], coords_min[0, 1], 'ko')
            axs[0].plot(coords_min[1, 0], coords_min[1, 1], 'ko')
        axs[0].set_aspect('equal')
        axs[0].set_title("x vs y", fontsize=25)

        axs[1].plot(i.x, i.z, ',')
        if coords_min.shape == (3,):
            axs[1].plot(coords_min[0], coords_min[2], 'ko')
        if coords_min.shape == (2, 3):
            axs[1].plot(coords_min[0, 0], coords_min[0, 2], 'ko')
            axs[1].plot(coords_min[1, 0], coords_min[1, 2], 'ko')
        axs[1].set_aspect('equal')
        axs[1].set_title("x vs z", fontsize=25)

        axs[2].plot(i.y, i.z, ',')
        if coords_min.shape == (3,):
            axs[2].plot(coords_min[1], coords_min[2], 'ko')
        if coords_min.shape == (2, 3):
            axs[2].plot(coords_min[0, 1], coords_min[0, 2], 'ko')
            axs[2].plot(coords_min[1, 1], coords_min[1, 2], 'ko')
        axs[2].set_aspect('equal')
        axs[2].set_title("y vs z", fontsize=25)
        plt.show()
