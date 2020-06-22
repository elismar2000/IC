import astropy.table


def read_table(t_path, snapshot):
    '''
    Read mass of recently formed stars
    and their x, y, z positions for each snapshot
    '''
    t = astropy.table.Table.read(t_path, snapshot)

    gas = t['P_TYPE'] == 0
    new_born = t['M_GAS'][gas] < 0.05 * t['MASS'][gas]

    x = t['X'][gas][new_born]
    y = t['Y'][gas][new_born]
    z = t['Z'][gas][new_born]

    mass = 0.95 * t['MASS'][gas][new_born]

    return x, y, z, mass

#===============================

t_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/Tabelas_GalMer/tables_arp245_orbit1'
snapshot = 71

x, y, z, mass = read_table(t_path, snapshot)
