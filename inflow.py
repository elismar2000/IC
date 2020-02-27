from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt
from potential_3d import Potential
import pickle
import os

class Inflow:
    def __init__(self, t_path: str, snapshot: int):
        self.t_path = t_path
        self.snapshot = snapshot
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.m = np.array([])
        self.gm = np.array([])

    def read_table(self, p_type: int = None, gal: int = None):
        t = Table.read(self.t_path, self.snapshot)

        if p_type == None:
            gal_mask = t['GAL'] == gal
            x = t['X'][gal_mask]
            y = t['Y'][gal_mask]
            z = t['Z'][gal_mask]
            m = t['MASS'][gal_mask]
            gm = t['M_GAS'][gal_mask]

        if gal == None:
            p_type_mask = t['P_TYPE'] == p_type
            x = t['X'][p_type_mask]
            y = t['Y'][p_type_mask]
            z = t['Z'][p_type_mask]
            m = t['MASS'][p_type_mask]
            gm = t['M_GAS'][p_type_mask]

        if (gal == None) & (p_type == None):
            x = t['X']
            y = t['Y']
            z = t['Z']
            m = t['MASS']
            gm = t['M_GAS']

        if (gal != None) & (p_type != None):
            p_type_mask = t['P_TYPE'] == p_type
            gal_mask = t['GAL'][p_type_mask] == gal
            x = t['X'][p_type_mask][gal_mask]
            y = t['Y'][p_type_mask][gal_mask]
            z = t['Z'][p_type_mask][gal_mask]
            m = t['MASS'][p_type_mask][gal_mask]
            gm = t['M_GAS'][p_type_mask][gal_mask]

        self.x = x
        self.y = y
        self.z = z
        self.m = m
        self.gm = gm

    def select(self, orbit, which_min, min_radius: float = 1.0):
        with open('orbit' + str(orbit) + '_mins.txt', 'rb') as f:
            table = pickle.load(f)

        if which_min == 1:
            center = table[self.snapshot - 1, 1]
        if which_min == 2:
            center = table[self.snapshot - 1, 2]
        if which_min == 3:
            center = table[self.snapshot - 1, 3]

        r = np.array([self.x, self.y, self.z])
        dist = np.array([np.linalg.norm(r[:, i] - center) for i in range(len(r[0]))])
        mask = dist <= min_radius

        self.x = self.x[mask]
        self.y = self.y[mask]
        self.z = self.z[mask]
        self.m = self.m[mask]
        self.gm = self.gm[mask]

    def mass(self, p_type: int, gal: int = None, min_radius: float = 1.0, orbit: int = 1, which_min: int = 1):
        #p_type = 0 --> gas
        #p_type = 1 --> estrelas
        #p_type = 2 --> matéria escura
        if gal == None:
            if p_type == 0:
                self.read_table(0)
                self.select(orbit, which_min, min_radius)
                m = np.sum(self.gm)

            if p_type == 1:
                self.read_table(0)
                self.select(orbit, which_min, min_radius)
                m_gas = self.gm
                m_hybrid = self.m

                self.read_table(1)
                self.select(orbit, which_min, min_radius)
                m_star = self.m
                m_formed_star = m_hybrid - m_gas
                m = np.sum(m_star) + np.sum(m_formed_star)

            if p_type == 2:
                self.read_table(2)
                self.select(orbit, which_min, min_radius)
                m = np.sum(self.m)

        if gal != None:
            if p_type == 0:
                self.read_table(0, gal)
                self.select(orbit, which_min, min_radius)
                m = np.sum(self.gm)

            if p_type == 1:
                self.read_table(0, gal)
                self.select(orbit, which_min, min_radius)
                m_gas = self.gm
                m_hybrid = self.m

                self.read_table(1, gal)
                self.select(orbit, which_min, min_radius)
                m_star = self.m
                m_formed_star = m_hybrid - m_gas
                m = np.sum(m_star) + np.sum(m_formed_star)

            if p_type == 2:
                self.read_table(2, gal)
                self.select(orbit, which_min, min_radius)
                m = np.sum(self.m)

        return m

    def gas_particles(self, min_radius: float = 1.0, orbit: int = 1, which_min: int = 1):
        self.read_table(p_type=0)
        self.select(orbit, which_min, min_radius)
        return np.size(self.x)



if __name__ == '__main__':
    t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp142_v2'
    snapshot = 10
    i = Inflow(t_path, snapshot)
    i.read_table()
    plt.plot(i.x, i.y, ',')
    plt.show()
    # mode = 'part_num'
    #
    # for orbit in [1, 5, 9]:
    #     if orbit == 1:
    #         merge = 17
    #     if orbit == 5:
    #         merge = 29
    #     if orbit == 9:
    #         merge = 71
    #
    #     print('orbit = ', orbit)
    #     print('merge = ', merge)
    #
    #     t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit' + str(orbit)
    #
    #     for min_radius in [1.0, 0.5, 0.1, 0.01]:
    #         os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow')
    #         print('min_radius = ', min_radius)
    #
    #         if mode == 'part_num':
    #             n1 = np.zeros(merge)
    #             n2 = np.zeros(merge)
    #             n = np.zeros(71 - merge)
    #
    #             for snapshot in range(1, 71):
    #                 i = Inflow(t_path, snapshot)
    #                 print('snapshot: ', snapshot)
    #                 if snapshot <= merge:
    #                     n1[snapshot - 1] = i.gas_particles(min_radius=min_radius, orbit=orbit, which_min=1)
    #                     n2[snapshot - 1] = i.gas_particles(min_radius=min_radius, orbit=orbit, which_min=2)
    #                 if snapshot > merge:
    #                     n[snapshot - merge - 1] = i.gas_particles(min_radius=min_radius, orbit=orbit, which_min=3)
    #
    #             dir = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs'
    #             os.chdir(dir)
    #             np.save('n1.orbit' + str(orbit) + '.min_radius' + str(min_radius), n1)
    #             np.save('n2.orbit' + str(orbit) + '.min_radius' + str(min_radius), n2)
    #             np.save('n.orbit' + str(orbit) + '.min_radius' + str(min_radius), n)
    #
    #
    #         if mode == 'inflow':
    #             m_gm1 = np.zeros(merge)
    #             m_gm2 = np.zeros(merge)
    #             m_gm = np.zeros(71 - merge)
    #
    #             for snapshot in range(1, 72):
    #                 i = Inflow(t_path, snapshot)
    #                 print('snapshot: ', snapshot)
    #                 if snapshot <= merge:
    #                     m_gm1[snapshot - 1] = i.mass(p_type=0, orbit=orbit, which_min=1, min_radius=min_radius)
    #                     m_gm2[snapshot - 1] = i.mass(p_type=0, orbit=orbit, which_min=2, min_radius=min_radius)
    #                 if snapshot > merge:
    #                     m_gm[snapshot - merge - 1] = i.mass(p_type=0, orbit=orbit, which_min=3, min_radius=min_radius)
    #
    #             dir = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs'
    #             os.chdir(dir)
    #             np.save('m_gm1.orbit' + str(orbit) + '.min_radius' + str(min_radius), m_gm1)
    #             np.save('m_gm2.orbit' + str(orbit) + '.min_radius' + str(min_radius), m_gm2)
    #             np.save('m_gm.orbit' + str(orbit) + '.min_radius' + str(min_radius), m_gm)
    #
    # os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow')
