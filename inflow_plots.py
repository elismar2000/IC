import numpy as np
import matplotlib.pyplot as plt
import os

cmap = plt.get_cmap('plasma')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 3, figsize=(40, 40), sharey='row', sharex=True)

f = 2.25e+9 #Fator da unidade do GalMer

from matplotlib import rc
import matplotlib.font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('text', usetex=True)

l = -1
for i in [1, 5, 9]:
    l += 1
    k = -1
    for j in [1.0, 0.5, 0.1, 0.01]:
        k += 1
        #v = 4/3 * np.pi * j**3

        if (i == 1) | (i == 5):
            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs')
            gm1 = np.load('m_gm1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm2 = np.load('m_gm2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm = np.load('m_gm.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs')
            n1 = np.load('n1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            n2 = np.load('n2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            n = np.load('n.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            flux1 = np.diff(gm1) / 50e+6
            flux2 = np.diff(gm2) / 50e+6
            flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

            err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
            err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
            err = (np.sqrt(n) * gm / n) / 50e+6

        if (i == 9):
            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs')
            gm1 = np.load('m_gm1.orbit' + str(i) + '.min_radius' + str(j) + '.npy')
            gm2 = np.load('m_gm2.orbit' + str(i) + '.min_radius' + str(j) + '.npy')

            os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs')
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
                axs[k, l].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5, label=r'Galaxy 1')
                axs[k, l].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5, label=r'Galaxy 2')
                axs[k, l].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5, label=r'Post merger galaxy')

            else:
                axs[k, l].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5)
                axs[k, l].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5)
                axs[k, l].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5)

            axs[k, l].fill_between(time1, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k, l].fill_between(time1, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
            axs[k, l].fill_between(time2, f * (flux - err), f * (flux + err), alpha=0.3, color=colors[5])
            axs[k, l].axvline(350e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)

            axs[k, l].set_ylabel(r'$\frac{M\odot}{yr}$', fontsize=25)

        if i == 5:
            time1 = np.arange(1, 29, 1)*0.05 #Gyr
            time2 = np.arange(29, 71, 1)*0.05 #Gyr

            axs[k, l].plot(time1, (flux1 * f), '--', color=colors[3], linewidth=1.5)
            axs[k, l].plot(time1, (flux2 * f), '--', color=colors[1], linewidth=1.5)
            axs[k, l].plot(time2, (flux * f), '--', color=colors[5], linewidth=1.5)
            axs[k, l].axvline(400e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k, l].axvline(1250e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k, l].fill_between(time1, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k, l].fill_between(time1, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
            axs[k, l].fill_between(time2, f * (flux - err), f * (flux + err), alpha=0.3, color=colors[5])

        if i == 9:
            time = np.arange(1, 71, 1)*0.05 #Gyr

            axs[k, l].plot(time, flux1 * f, '--', color=colors[3], linewidth=1.5)
            axs[k, l].plot(time, flux2 * f, '--', color=colors[1], linewidth=1.5)
            axs[k, l].axvline(450e-3, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
            axs[k, l].fill_between(time, f * (flux1 - err1[1:]), f * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
            axs[k, l].fill_between(time, f * (flux2 - err2[1:]), f * (flux2 + err2[1:]), alpha=0.3, color=colors[1])

            axs[k, l].yaxis.set_label_position('right')
            axis_cmap = plt.get_cmap('plasma')
            axis_colors = [cmap(i) for i in np.linspace(0, 1, 7)]
            if j == 1.0:
                axs[k, l].yaxis.label.set_color(axis_colors[0])
                axs[k, l].set_ylabel(r'$r = 1.0 kpc$', fontsize=20)
            if j == 0.5:
                axs[k, l].yaxis.label.set_color(axis_colors[3])
                axs[k, l].set_ylabel(r'$r = 0.5 kpc$', fontsize=20)
            if j == 0.1:
                axs[k, l].yaxis.label.set_color(axis_colors[1])
                axs[k, l].set_ylabel(r'$r = 0.1 kpc$', fontsize=20)
            if j == 0.01:
                axs[k, l].yaxis.label.set_color(axis_colors[2])
                axs[k, l].set_ylabel(r'$r = 0.01 kpc$', fontsize=20)

        if j == 1.0:
            axs[k, l].xaxis.set_label_position('top')
            axis_cmap = plt.get_cmap('rainbow')
            axis_colors = [cmap(i) for i in np.linspace(0, 1, 7)]
            if i == 1:
                axs[k, l].xaxis.label.set_color(axis_colors[3])
                axs[k, l].set_xlabel(r'$Orbit\ type\ 1$', fontsize=30)
            if i == 5:
                axs[k, l].xaxis.label.set_color(axis_colors[4])
                axs[k, l].set_xlabel(r'$Orbit\ type\ 5$', fontsize=30)
            if i == 9:
                axs[k, l].xaxis.label.set_color(axis_colors[1])
                axs[k, l].set_xlabel(r'$Orbit\ type\ 9$', fontsize=30)

        axs[k, l].legend(loc=0, fontsize=10, fancybox=False, framealpha=0.0)
        axs[k, l].tick_params(axis='both', which='both', tickdir='in', labelsize=25,
            top=True, right=True, width=1.5)
        axs[k, l].tick_params(which='minor', length=3)
        axs[k, l].tick_params(which='major', length=6)
        axs[k, l].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axs[k, l].minorticks_on()
        axs[k, l].spines['bottom'].set_linewidth(1.5)
        axs[k, l].spines['left'].set_linewidth(1.5)
        axs[k, l].spines['top'].set_linewidth(1.5)
        axs[k, l].spines['right'].set_linewidth(1.5)

axs[3, 1].set_xlabel(r'$Gyr$', fontsize=30)

plt.subplots_adjust(wspace=0, hspace=0)

os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow')

plt.savefig('gas_flux.pdf', dpi=0.1)
plt.show()
