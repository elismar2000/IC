import matplotlib.pyplot as plt
import numpy as np

from tables import *

#===========================
#Reading tables
#===========================

sfr_path_gadget = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/sfr_outputs_gadget/galmer-like_sim_sfr.h5'
sfr_path_galmer = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/sfr_outputs_galmer/orbit1_galmer_sim_sfr.h5'

sfr_h5file_gadget = open_file(sfr_path_gadget, 'r')
sfr_table_gadget = sfr_h5file_gadget.root.sfr.readout

sfr_h5file_galmer = open_file(sfr_path_galmer, 'r')
sfr_table_galmer = sfr_h5file_galmer.root.sfr.readout

#===========================
#Plotting
#===========================

cmap = plt.get_cmap('twilight')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharey=False, sharex=True, dpi=60)

#============================
#unit factors
#============================

age_factor_gg = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor_gg = 1e+10

age_factor_gm = 50e+6 #Time between each snapshot of the simulation in years
mass_factor_gm = 2.25e+9

time = np.arange(1, 71, 1) * 0.05 #Time in Gyr

#============================
#Radius = 1.2 Kpc
#============================
#Gadget
#============================

gm1 = np.array([ j['gal1_r10'] for j in sfr_table_gadget.iterrows()])
gm2 = np.array([ j['gal2_r10'] for j in sfr_table_gadget.iterrows()])
gm = np.array([ j['gal_r10'] for j in sfr_table_gadget.iterrows()])

axs[0].plot(time, np.log10(gm1 * mass_factor_gg), '-', color=colors[1], linewidth=3.0)
axs[0].plot(time, np.log10(gm2 * mass_factor_gg), '-', color=colors[3], linewidth=3.0)
axs[0].plot(time, np.log10(gm * mass_factor_gg), '-', color=colors[5], linewidth=3.0)

#============================
#GalMer
#============================

gm1 = np.array([ j['gal1_r10'] for j in sfr_table_galmer.iterrows()])
gm2 = np.array([ j['gal2_r10'] for j in sfr_table_galmer.iterrows()])
gm = np.array([ j['gal_r10'] for j in sfr_table_galmer.iterrows()])

axs[0].plot(time, np.log10(gm1 * mass_factor_gm), '--', alpha = 0.4, color=colors[1], linewidth=3.0)
axs[0].plot(time, np.log10(gm2 * mass_factor_gm), '--', alpha = 0.4, color=colors[3], linewidth=3.0)
axs[0].plot(time, np.log10(gm * mass_factor_gm), '--', alpha = 0.4, color=colors[5], linewidth=3.0)

#============================

axs[0].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=False, right=False, width=1.5)
axs[0].tick_params(which='minor', length=3)
axs[0].tick_params(which='major', length=6)
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[0].minorticks_on()

axs[0].text(0.2, 0.62*1e1, 'r = 1.0 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.5 Kpc
#============================
#Gadget
#============================

gm1 = np.array([ j['gal1_r05'] for j in sfr_table_gadget.iterrows()])
gm2 = np.array([ j['gal2_r05'] for j in sfr_table_gadget.iterrows()])
gm = np.array([ j['gal_r05'] for j in sfr_table_gadget.iterrows()])

axs[1].plot(time, np.log10(gm1 * mass_factor_gg), '-', color=colors[1], linewidth=3.0)
axs[1].plot(time, np.log10(gm2 * mass_factor_gg), '-', color=colors[3], linewidth=3.0)
axs[1].plot(time, np.log10(gm * mass_factor_gg), '-', color=colors[5], linewidth=3.0)

#============================
#GalMer
#============================

gm1 = np.array([ j['gal1_r05'] for j in sfr_table_galmer.iterrows()])
gm2 = np.array([ j['gal2_r05'] for j in sfr_table_galmer.iterrows()])
gm = np.array([ j['gal_r05'] for j in sfr_table_galmer.iterrows()])

axs[1].plot(time, np.log10(gm1 * mass_factor_gm), '--', alpha = 0.4, color=colors[1], linewidth=3.0)
axs[1].plot(time, np.log10(gm2 * mass_factor_gm), '--', alpha = 0.4, color=colors[3], linewidth=3.0)
axs[1].plot(time, np.log10(gm * mass_factor_gm), '--', alpha = 0.4, color=colors[5], linewidth=3.0)

#============================

axs[1].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=False, right=False, width=1.5)
axs[1].tick_params(which='minor', length=3)
axs[1].tick_params(which='major', length=6)
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[1].minorticks_on()

axs[1].text(0.2, 0.59*1e1, 'r = 0.5 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.1 Kpc
#============================
#Gadget
#============================

gm1 = np.array([ j['gal1_r01'] for j in sfr_table_gadget.iterrows()])
gm2 = np.array([ j['gal2_r01'] for j in sfr_table_gadget.iterrows()])
gm = np.array([ j['gal_r01'] for j in sfr_table_gadget.iterrows()])

axs[2].plot(time, np.log10(gm1 * mass_factor_gg), '-', color=colors[1], linewidth=3.0)
axs[2].plot(time, np.log10(gm2 * mass_factor_gg), '-', color=colors[3], linewidth=3.0)
axs[2].plot(time, np.log10(gm * mass_factor_gg), '-', color=colors[5], linewidth=3.0)

#============================
#GalMer
#============================

gm1 = np.array([ j['gal1_r01'] for j in sfr_table_galmer.iterrows()])
gm2 = np.array([ j['gal2_r01'] for j in sfr_table_galmer.iterrows()])
gm = np.array([ j['gal_r01'] for j in sfr_table_galmer.iterrows()])

axs[2].plot(time, np.log10(gm1 * mass_factor_gm), '--', alpha = 0.4, color=colors[1], linewidth=3.0)
axs[2].plot(time, np.log10(gm2 * mass_factor_gm), '--', alpha = 0.4, color=colors[3], linewidth=3.0)
axs[2].plot(time, np.log10(gm * mass_factor_gm), '--', alpha = 0.4, color=colors[5], linewidth=3.0)

#============================

axs[2].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=False, right=False, width=1.5)
axs[2].tick_params(which='minor', length=3)
axs[2].tick_params(which='major', length=6)
axs[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[2].minorticks_on()

axs[2].text(0.2, 0.58*1e1, 'r = 0.1 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.01 Kpc
#============================
#Gadget
#============================

gm1 = np.array([ j['gal1_r001'] for j in sfr_table_gadget.iterrows()])
gm2 = np.array([ j['gal2_r001'] for j in sfr_table_gadget.iterrows()])
gm = np.array([ j['gal_r001'] for j in sfr_table_gadget.iterrows()])

axs[3].plot(time, np.log10(gm1 * mass_factor_gg), '-', color=colors[1], linewidth=3.0)
axs[3].plot(time, np.log10(gm2 * mass_factor_gg), '-', color=colors[3], linewidth=3.0)
axs[3].plot(time, np.log10(gm * mass_factor_gg), '-', color=colors[5], linewidth=3.0)

#============================
#GalMer
#============================

gm1 = np.array([ j['gal1_r001'] for j in sfr_table_galmer.iterrows()])
gm2 = np.array([ j['gal2_r001'] for j in sfr_table_galmer.iterrows()])
gm = np.array([ j['gal_r001'] for j in sfr_table_galmer.iterrows()])

axs[3].plot(time, np.log10(gm1 * mass_factor_gm), '--', alpha = 0.4, color=colors[1], linewidth=3.0)
axs[3].plot(time, np.log10(gm2 * mass_factor_gm), '--', alpha = 0.4, color=colors[3], linewidth=3.0)
axs[3].plot(time, np.log10(gm * mass_factor_gm), '--', alpha = 0.4, color=colors[5], linewidth=3.0)

#============================

axs[3].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=False, right=False, width=1.5)
axs[3].tick_params(which='minor', length=3)
axs[3].tick_params(which='major', length=6)
axs[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[3].minorticks_on()

axs[3].text(0.2, 0.19*1e1, 'r = 0.01 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================

axs[3].plot([1, 1],[1, 1], '-', color='k', label='Gadget simulation')
axs[3].plot([1, 1],[1, 1], '--', color='k', label='Galmer simulation')
axs[3].plot([1, 1],[1, 1], 'o', color=colors[1], label='Galaxy 1')
axs[3].plot([1, 1],[1, 1], 'o', color=colors[3], label='Galaxy 2')
axs[3].plot([1, 1],[1, 1], 'o', color=colors[5], label='Coalesced galaxy')
axs[3].plot([1, 1],[1, 1], 'o', color='white')
axs[3].legend(loc=0, fontsize=10, fancybox=True, framealpha=0.4, shadow=True)

fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.ylabel(r'$log(SFR)\ [(\frac{M\odot}{yr})]$', fontsize=20)
plt.xlabel(r'$Time\ [Gyr]$', fontsize=20)

plt.subplots_adjust(wspace=0, hspace=0)

import os
os.chdir('/home/elismar/Documentos/Fisica/IC/figuras/RAS Poster')
plt.savefig('sfr_plot.png')
#plt.show()

os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/codes')
