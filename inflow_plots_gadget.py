import matplotlib.pyplot as plt
import numpy as np

from tables import *

#===========================
#Reading tables
#===========================

inflow_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_gadget/galmer-like_sim_inflow.h5'
part_num_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_gadget/galmer-like_sim_part_num.h5'

inflow_h5file = open_file(inflow_path, 'r')
inflow_table = inflow_h5file.root.inflow.readout

part_num_h5file = open_file(part_num_path, 'r')
part_num_table = part_num_h5file.root.part_num.readout

#===========================
#Plotting
#===========================

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharey=False, sharex=True)

#============================
#Radius = 1.0 Kpc
#============================

gm1_r10 = np.array([ j['gal1_r10'] for j in inflow_table.iterrows()])
gm2_r10 = np.array([ j['gal2_r10'] for j in inflow_table.iterrows()])
gm_r10 = np.array([ j['gal_r10'] for j in inflow_table.iterrows()])

n1_r10 = np.array([ j['gal1_r10'] for j in part_num_table.iterrows()])
n2_r10 = np.array([ j['gal2_r10'] for j in part_num_table.iterrows()])
n_r10 = np.array([ j['gal_r10'] for j in part_num_table.iterrows()])

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

flux_gm1_r10 = np.diff(gm1_r10) / age_factor
flux_gm2_r10 = np.diff(gm2_r10) / age_factor
flux_gm_r10 = np.diff(gm_r10) / age_factor

err1 = (np.sqrt(n1_r10) * gm1_r10 / n1_r10) / age_factor
err2 = (np.sqrt(n2_r10) * gm2_r10 / n2_r10) / age_factor
err = (np.sqrt(n_r10) * gm_r10 / n_r10) / age_factor

time = np.arange(1, 71, 1) * 0.05 * 0.98 #Time in Gyr

axs[0].plot(time, flux_gm1_r10 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[0].plot(time, flux_gm2_r10 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[0].plot(time, flux_gm_r10 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[0].fill_between(time, mass_factor * (flux_gm1_r10 - err1[1:]), mass_factor * (flux_gm1_r10 + err1[1:]), alpha=0.3, color='navy')
axs[0].fill_between(time, mass_factor * (flux_gm2_r10 - err2[1:]), mass_factor * (flux_gm2_r10 + err2[1:]), alpha=0.3, color='darkorange')
axs[0].fill_between(time, mass_factor * (flux_gm_r10 - err[1:]), mass_factor * (flux_gm_r10 + err[1:]), alpha=0.3, color='m')

axs[0].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[0].tick_params(which='minor', length=3)
axs[0].tick_params(which='major', length=6)
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[0].minorticks_on()

axs[0].set_ylabel(r'$Flux\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.5 Kpc
#============================

gm1_r05 = np.array([ j['gal1_r05'] for j in inflow_table.iterrows()])
gm2_r05 = np.array([ j['gal2_r05'] for j in inflow_table.iterrows()])
gm_r05 = np.array([ j['gal_r05'] for j in inflow_table.iterrows()])

n1_r05 = np.array([ j['gal1_r05'] for j in part_num_table.iterrows()])
n2_r05 = np.array([ j['gal2_r05'] for j in part_num_table.iterrows()])
n_r05 = np.array([ j['gal_r05'] for j in part_num_table.iterrows()])

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

flux_gm1_r05 = np.diff(gm1_r05) / age_factor
flux_gm2_r05 = np.diff(gm2_r05) / age_factor
flux_gm_r05 = np.diff(gm_r05) / age_factor

err1 = (np.sqrt(n1_r05) * gm1_r05 / n1_r05) / age_factor
err2 = (np.sqrt(n2_r05) * gm2_r05 / n2_r05) / age_factor
err = (np.sqrt(n_r05) * gm_r05 / n_r05) / age_factor

time = np.arange(1, 71, 1) * 0.05 * 0.98 #Time in Gyr

axs[1].plot(time, flux_gm1_r05 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[1].plot(time, flux_gm2_r05 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[1].plot(time, flux_gm_r05 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[1].fill_between(time, mass_factor * (flux_gm1_r05 - err1[1:]), mass_factor * (flux_gm1_r05 + err1[1:]), alpha=0.3, color='navy')
axs[1].fill_between(time, mass_factor * (flux_gm2_r05 - err2[1:]), mass_factor * (flux_gm2_r05 + err2[1:]), alpha=0.3, color='darkorange')
axs[1].fill_between(time, mass_factor * (flux_gm_r05 - err[1:]), mass_factor * (flux_gm_r05 + err[1:]), alpha=0.3, color='m')

axs[1].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[1].tick_params(which='minor', length=3)
axs[1].tick_params(which='major', length=6)
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[1].minorticks_on()

axs[1].set_ylabel(r'$Flux\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.1 Kpc
#============================

gm1_r01 = np.array([ j['gal1_r01'] for j in inflow_table.iterrows()])
gm2_r01 = np.array([ j['gal2_r01'] for j in inflow_table.iterrows()])
gm_r01 = np.array([ j['gal_r01'] for j in inflow_table.iterrows()])

n1_r01 = np.array([ j['gal1_r01'] for j in part_num_table.iterrows()])
n2_r01 = np.array([ j['gal2_r01'] for j in part_num_table.iterrows()])
n_r01 = np.array([ j['gal_r01'] for j in part_num_table.iterrows()])

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

flux_gm1_r01 = np.diff(gm1_r01) / age_factor
flux_gm2_r01 = np.diff(gm2_r01) / age_factor
flux_gm_r01 = np.diff(gm_r01) / age_factor

err1 = (np.sqrt(n1_r01) * gm1_r01 / n1_r01) / age_factor
err2 = (np.sqrt(n2_r01) * gm2_r01 / n2_r01) / age_factor
err = (np.sqrt(n_r01) * gm_r01 / n_r01) / age_factor

time = np.arange(1, 71, 1) * 0.05 * 0.98 #Time in Gyr

axs[2].plot(time, flux_gm1_r01 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[2].plot(time, flux_gm2_r01 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[2].plot(time, flux_gm_r01 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[2].fill_between(time, mass_factor * (flux_gm1_r01 - err1[1:]), mass_factor * (flux_gm1_r01 + err1[1:]), alpha=0.3, color='navy')
axs[2].fill_between(time, mass_factor * (flux_gm2_r01 - err2[1:]), mass_factor * (flux_gm2_r01 + err2[1:]), alpha=0.3, color='darkorange')
axs[2].fill_between(time, mass_factor * (flux_gm_r01 - err[1:]), mass_factor * (flux_gm_r01 + err[1:]), alpha=0.3, color='m')

axs[2].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[2].tick_params(which='minor', length=3)
axs[2].tick_params(which='major', length=6)
axs[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[2].minorticks_on()

axs[2].set_ylabel(r'$Flux\ [\frac{M\odot}{yr}]$', fontsize=12)

#============================
#Radius = 0.01 Kpc
#============================

gm1_r001 = np.array([ j['gal1_r001'] for j in inflow_table.iterrows()])
gm2_r001 = np.array([ j['gal2_r001'] for j in inflow_table.iterrows()])
gm_r001 = np.array([ j['gal_r001'] for j in inflow_table.iterrows()])

n1_r001 = np.array([ j['gal1_r001'] for j in part_num_table.iterrows()])
n2_r001 = np.array([ j['gal2_r001'] for j in part_num_table.iterrows()])
n_r001 = np.array([ j['gal_r001'] for j in part_num_table.iterrows()])

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

flux_gm1_r001 = np.diff(gm1_r001) / age_factor
flux_gm2_r001 = np.diff(gm2_r001) / age_factor
flux_gm_r001 = np.diff(gm_r001) / age_factor

err1 = (np.sqrt(n1_r001) * gm1_r001 / n1_r001) / age_factor
err2 = (np.sqrt(n2_r001) * gm2_r001 / n2_r001) / age_factor
err = (np.sqrt(n_r001) * gm_r001 / n_r001) / age_factor

time = np.arange(1, 71, 1) * 0.05 * 0.98 #Time in Gyr

axs[3].plot(time, flux_gm1_r001 * mass_factor, '-', color='navy', linewidth=1.5, label=r'Galaxy 1')
axs[3].plot(time, flux_gm2_r001 * mass_factor, '-', color='darkorange', linewidth=1.5, label=r'Galaxy 2')
axs[3].plot(time, flux_gm_r001 * mass_factor, '-', color='m', linewidth=1.5, label=r'Galaxy')

axs[3].fill_between(time, mass_factor * (flux_gm1_r001 - err1[1:]), mass_factor * (flux_gm1_r001 + err1[1:]), alpha=0.3, color='navy')
axs[3].fill_between(time, mass_factor * (flux_gm2_r001 - err2[1:]), mass_factor * (flux_gm2_r001 + err2[1:]), alpha=0.3, color='darkorange')
axs[3].fill_between(time, mass_factor * (flux_gm_r001 - err[1:]), mass_factor * (flux_gm_r001 + err[1:]), alpha=0.3, color='m')

axs[3].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[3].tick_params(which='minor', length=3)
axs[3].tick_params(which='major', length=6)
axs[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[3].minorticks_on()

axs[3].set_ylabel(r'$Flux\ [\frac{M\odot}{yr}]$', fontsize=12)
axs[3].set_xlabel(r'$Time\ [Gyr]$', fontsize=15)

#============================

plt.legend(loc=0, fontsize=10, fancybox=True, framealpha=0.4, shadow=True)

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
