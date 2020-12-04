import matplotlib.pyplot as plt
import numpy as np

from tables import *

#===========================
#Reading tables
#===========================

sfr_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/sfr_outputs_galmer/orbit1_galmer_sim_sfr.h5'

sfr_h5file = open_file(sfr_path, 'r')
sfr_table = sfr_h5file.root.sfr.readout

#===========================
#Plotting
#===========================

cmap = plt.get_cmap('twilight')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharey=False, sharex=True)

#============================
#Radius = 1.0 Kpc
#============================

gm1_r10 = np.array([ j['gal1_r10'] for j in sfr_table.iterrows()])
gm2_r10 = np.array([ j['gal2_r10'] for j in sfr_table.iterrows()])
gm_r10 = np.array([ j['gal_r10'] for j in sfr_table.iterrows()])

age_factor = 50e+6 #Time between each snapshot of the simulation in years
mass_factor = 2.25e+9

time = np.arange(1, 71, 1) * 0.05 #Time in Gyr

axs[0].plot(time, gm1_r10 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[0].plot(time, gm2_r10 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[0].plot(time, gm_r10 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[0].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[0].tick_params(which='minor', length=3)
axs[0].tick_params(which='major', length=6)
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[0].minorticks_on()

axs[0].set_ylabel(r'$SFR\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.5 Kpc
#============================

gm1_r05 = np.array([ j['gal1_r05'] for j in sfr_table.iterrows()])
gm2_r05 = np.array([ j['gal2_r05'] for j in sfr_table.iterrows()])
gm_r05 = np.array([ j['gal_r05'] for j in sfr_table.iterrows()])

age_factor = 50e+6 #Time between each snapshot of the simulation in years
mass_factor = 2.25e+9

time = np.arange(1, 71, 1) * 0.05 #Time in Gyr

axs[1].plot(time, gm1_r05 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[1].plot(time, gm2_r05 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[1].plot(time, gm_r05 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[1].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[1].tick_params(which='minor', length=3)
axs[1].tick_params(which='major', length=6)
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[1].minorticks_on()

axs[1].set_ylabel(r'$SFR\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.1 Kpc
#============================

gm1_r01 = np.array([ j['gal1_r01'] for j in sfr_table.iterrows()])
gm2_r01 = np.array([ j['gal2_r01'] for j in sfr_table.iterrows()])
gm_r01 = np.array([ j['gal_r01'] for j in sfr_table.iterrows()])

age_factor = 50e+6 #Time between each snapshot of the simulation in years
mass_factor = 2.25e+9

time = np.arange(1, 71, 1) * 0.05 #Time in Gyr

axs[2].plot(time, gm1_r01 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[2].plot(time, gm2_r01 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[2].plot(time, gm_r01 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[2].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[2].tick_params(which='minor', length=3)
axs[2].tick_params(which='major', length=6)
axs[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[2].minorticks_on()

axs[2].set_ylabel(r'$SFR\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.01 Kpc
#============================

gm1_r001 = np.array([ j['gal1_r001'] for j in sfr_table.iterrows()])
gm2_r001 = np.array([ j['gal2_r001'] for j in sfr_table.iterrows()])
gm_r001 = np.array([ j['gal_r001'] for j in sfr_table.iterrows()])

age_factor = 50e+6 #Time between each snapshot of the simulation in years
mass_factor = 2.25e+9

time = np.arange(1, 71, 1) * 0.05 #Time in Gyr

axs[3].plot(time, gm1_r001 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[3].plot(time, gm2_r001 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[3].plot(time, gm_r001 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[3].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[3].tick_params(which='minor', length=3)
axs[3].tick_params(which='major', length=6)
axs[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[3].minorticks_on()

axs[3].set_ylabel(r'$SFR\ [\frac{M\odot}{yr}]$', fontsize=12)
axs[3].set_xlabel(r'$Time\ [Gyr]$', fontsize=15)

#============================

plt.legend(loc=0, fontsize=10, fancybox=True, framealpha=0.4, shadow=True)

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
