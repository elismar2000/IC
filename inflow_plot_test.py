import numpy as np
import matplotlib.pyplot as plt
import os

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharey='row', sharex=True)

f = 2.25e+9 #Fator da unidade do GalMer

# from matplotlib import rc
# import matplotlib.font_manager
# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
# rc('text', usetex=True)

#l = -1
for i in [1]:
    #l += 1
    k = -1
    for j in [1.0, 0.5, 0.1, 0.01]:
        k += 1

        if (i == 1) | (i == 5):
            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer')
            gm1 = np.load('m_gm1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm2 = np.load('m_gm2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm = np.load('m_gm.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer')
            n1 = np.load('n1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            n2 = np.load('n2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            n = np.load('n.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            #50e+6 is a factor to let the flux in units of solar mass / year
            flux1 = np.diff(gm1) / 50e+6
            flux2 = np.diff(gm2) / 50e+6
            flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

            err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
            err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
            err = (np.sqrt(n) * gm / n) / 50e+6

        if (i == 9):
            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer')
            gm1 = np.load('m_gm1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm2 = np.load('m_gm2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer')
            n1 = np.load('n1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            n2 = np.load('n2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            flux1 = np.diff(gm1) / 50e+6
            flux2 = np.diff(gm2) / 50e+6

            err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
            err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6

        if i == 1:
            time1 = np.arange(1, 17, 1)*0.05 #Gyr
            time2 = np.arange(17, 71, 1)*0.05 #Gyr

            if j == 1.0:
                axs[k].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5, label=r'Galaxy 1')
                axs[k].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5, label=r'Galaxy 2')
                axs[k].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5, label=r'Post merger galaxy')

            else:
                axs[k].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5)
                axs[k].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5)
                axs[k].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5)

            axs[k].fill_between(time1, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k].fill_between(time1, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
            axs[k].fill_between(time2, f * (flux - err), f * (flux + err), alpha=0.3, color=colors[5])
            axs[k].axvline(350e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)

            axs[k].set_ylabel(r'$\frac{M\odot}{yr}$', fontsize=25)

        if i == 5:
            time1 = np.arange(1, 29, 1)*0.05 #Gyr
            time2 = np.arange(29, 71, 1)*0.05 #Gyr

            axs[k].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5)
            axs[k].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5)
            axs[k].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5)
            axs[k].axvline(400e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k].axvline(1250e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k].fill_between(time1, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k].fill_between(time1, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
            axs[k].fill_between(time2, f * (flux - err), f * (flux + err), alpha=0.3, color=colors[5])

        if i == 9:
            time = np.arange(1, 71, 1)*0.05 #Gyr

            axs[k].plot(time, flux1 * f, '--', color=colors[3], linewidth=1.5)
            axs[k].plot(time, flux2 * f, '--', color=colors[1], linewidth=1.5)
            axs[k].axvline(450e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k].fill_between(time, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k].fill_between(time, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])

            axs[k, l].yaxis.set_label_position('right')
            axis_cmap = plt.get_cmap('plasma')
            axis_colors = [cmap(i) for i in np.linspace(0, 1, 7)]
            if j == 1.0:
                axs[k].yaxis.label.set_color('black')
                axs[k].set_ylabel(r'$1.0\ kpc$', fontsize=20)
            if j == 0.5:
                axs[k].yaxis.label.set_color('black')
                axs[k].set_ylabel(r'$0.5\ kpc$', fontsize=20)
            if j == 0.1:
                axs[k].yaxis.label.set_color('black')
                axs[k].set_ylabel(r'$0.1\ kpc$', fontsize=20)
            if j == 0.01:
                axs[k].yaxis.label.set_color('black')
                axs[k].set_ylabel(r'$0.01\ kpc$', fontsize=20)

        if j == 1.0:
            axs[k].xaxis.set_label_position('top')
            axis_cmap = plt.get_cmap('rainbow')
            axis_colors = [cmap(i) for i in np.linspace(0, 1, 7)]
            if i == 1:
                axs[k].xaxis.label.set_color(axis_colors[3])
                axs[k].set_xlabel(r'$Orbit\ type\ 1$', fontsize=25)
            if i == 5:
                axs[k].xaxis.label.set_color(axis_colors[4])
                axs[k].set_xlabel(r'$Orbit\ type\ 5$', fontsize=25)
            if i == 9:
                axs[k].xaxis.label.set_color(axis_colors[1])
                axs[k].set_xlabel(r'$Orbit\ type\ 9$', fontsize=25)

        axs[k].legend(loc=0, fontsize=10, fancybox=False, framealpha=0.0)
        axs[k].tick_params(axis='both', which='both', tickdir='in', labelsize=25,
            top=True, right=True, width=1.5)
        axs[k].tick_params(which='minor', length=3)
        axs[k].tick_params(which='major', length=6)
        axs[k].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axs[k].minorticks_on()
        axs[k].spines['bottom'].set_linewidth(1.5)
        axs[k].spines['left'].set_linewidth(1.5)
        axs[k].spines['top'].set_linewidth(1.5)
        axs[k].spines['right'].set_linewidth(1.5)

axs[3].set_xlabel(r'$Gyr$', fontsize=23)

plt.subplots_adjust(wspace=0, hspace=0)

#===========================
#plot do gadget
#===========================

from tables import *

#===========================
#Reading tables
#===========================

sfr_path = '/home/elismar/Documentos/Fisica/IC/GalMer/inflow/sfr_outputs_gadget/galmer-like_sim_sfr.h5'

sfr_h5file = open_file(sfr_path, 'r')
sfr_table = sfr_h5file.root.sfr.readout

#===========================
#Plotting
#===========================

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

#============================
#Radius = 1.0 Kpc
#============================

gm1_r10 = np.array([ j['gal1_r10'] for j in sfr_table.iterrows()])
gm2_r10 = np.array([ j['gal2_r10'] for j in sfr_table.iterrows()])
gm_r10 = np.array([ j['gal_r10'] for j in sfr_table.iterrows()])

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

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

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

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

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

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

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

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

plt.show()
