import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons
import os
from astropy import units

for i in [1, 5, 9]:
    os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs')
    gm1 = np.load('m_gm1.orbit' + str(i) + '.min_radius0.01.npy') * 2.25e+9
    gm2 = np.load('m_gm2.orbit' + str(i) + '.min_radius0.01.npy') * 2.25e+9

    os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs')
    n1 = np.load('n1.orbit' + str(i) + '.min_radius0.01.npy')
    n2 = np.load('n2.orbit' + str(i) + '.min_radius0.01.npy')

    flux1 = np.diff(gm1)
    flux2 = np.diff(gm2)

    err1 = np.sqrt(n1) * gm1 / n1
    err2 = np.sqrt(n2) * gm2 / n2

    flux1 -= err1[1:]
    flux2 -= err2[1:]

    mdot1 = np.zeros_like(flux1)
    mdot2 = np.zeros_like(flux2)

    for j in range(len(mdot1)):
        if flux1[j] > 0:
            mdot1[j] = flux1[j] / 50e+6
        if flux1[j] <= 0:
            mdot1[j] = 0
        if flux2[j] > 0:
            mdot2[j] = flux2[j] / 50e+6
        if flux2[j] <= 0:
            mdot2[j] = 0

    eff = 0.1
    #Luminosidade vai ser em luminosidades solares.
    if i == 1:
        lum1_orbit1 = np.log10(mdot1 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)
        lum2_orbit1 = np.log10(mdot2 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)

    if i == 5:
        lum1_orbit5 =  np.log10(mdot1 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)
        lum2_orbit5 =  np.log10(mdot2 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)

    if i == 9:
        lum1_orbit9 =  np.log10(mdot1 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)
        lum2_orbit9 =  np.log10(mdot2 * cons.c**2 * eff * 2e+30 / (3.154e+7 * 3.829e+26) / 1e+9)

flux_2992_xray = 4.820000e-11 #erg / s cm²
d = 9.25703274e+25 #cm
lum_2992_xray = 4 * np.pi * d**2 * flux_2992_xray #Lsu
lum_2992_bol = 10**(0.0378 * (np.log10(lum_2992_xray))**2 - 2.03 * np.log10(lum_2992_xray) + 61.6)
lum_2992_to_plot = np.log10(lum_2992_bol / (3.846e+33 * 1e+9))

cmap = plt.get_cmap('gnuplot')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

axis_cmap = plt.get_cmap('rainbow')
axis_colors = [cmap(i) for i in np.linspace(0, 1, 7)]

from matplotlib import rc
import matplotlib.font_manager
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
rc('text', usetex=True)

time1 = np.arange(1, 17, 1)*0.05 #Gyr
time5 = np.arange(1, 29, 1)*0.05 #Gyr
time9 = np.arange(1, 71, 1)*0.05 #Gyr

fig, ax = plt.subplots(1, 3, figsize=(20, 20), sharey=True)
ax[0].plot(time1, lum1_orbit1, 'X', color=colors[3], label=r'Galaxy 1')
ax[0].plot(time1, lum2_orbit1, 'X', color=colors[1], label=r'Galaxy 2')
ax[0].axhline(lum_2992_to_plot, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
ax[0].tick_params(axis='both', tickdir='in', length=10, labelsize=20, width=1.5, top=True, right=True)
ax[0].spines['bottom'].set_linewidth(1.5)
ax[0].spines['left'].set_linewidth(1.5)
ax[0].spines['top'].set_linewidth(1.5)
ax[0].spines['right'].set_linewidth(1.5)
ax[0].tick_params(which='minor', tickdir='in', length=3, top=True, right=True)
ax[0].tick_params(which='major', length=6)
ax[0].minorticks_on()
ax[0].set_ylabel(r'$log_{10}(10^{9}L_{\odot})$', fontsize=25)
ax[0].xaxis.set_label_position('top')
ax[0].xaxis.label.set_color(axis_colors[3])
ax[0].set_xlabel(r'$Orbit\ type\ 1$', fontsize=35)
ax[0].legend(loc=4, fontsize=23, fancybox=True, framealpha=0.8, shadow=True)
ax[0].set_ylim(-0.7, 2.4)
ax[0].axvline(350e-3, color='b', linestyle='dashed', linewidth=2.0, alpha=0.6)

ax[1].plot(time5, lum1_orbit5, 'X', color=colors[3])
ax[1].plot(time5, lum2_orbit5, 'X', color=colors[1])
ax[1].axhline(lum_2992_to_plot, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
ax[1].tick_params(axis='both', tickdir='in', length=10, labelsize=20, width=1.5, top=True, right=True)
ax[1].spines['bottom'].set_linewidth(1.5)
ax[1].spines['left'].set_linewidth(1.5)
ax[1].spines['top'].set_linewidth(1.5)
ax[1].spines['right'].set_linewidth(1.5)
ax[1].tick_params(which='minor', tickdir='in', length=3, top=True, right=True)
ax[1].tick_params(which='major', length=6)
ax[1].minorticks_on()
ax[1].set_xlabel(r'$Gyr$', fontsize=25)
ax[1].set_title(r'$Orbit\ type\ 5$', fontsize=35, color=axis_colors[4])
ax[1].set_ylim(-0.7, 2.4)
ax[1].axvline(400e-3, color='b', linestyle='dashed', linewidth=2.0, alpha=0.6)
ax[1].axvline(1250e-3, color='b', linestyle='dashed', linewidth=2.0, alpha=0.6)

ax[2].plot(time9, lum1_orbit9, 'X', color=colors[3])
ax[2].plot(time9, lum2_orbit9, 'X', color=colors[1])
ax[2].axhline(lum_2992_to_plot, color='k', linestyle='dashed', linewidth=2.0, alpha=0.6)
ax[2].tick_params(axis='both', tickdir='in', length=10, labelsize=20, width=1.5, top=True, right=True)
ax[2].spines['bottom'].set_linewidth(1.5)
ax[2].spines['left'].set_linewidth(1.5)
ax[2].spines['top'].set_linewidth(1.5)
ax[2].spines['right'].set_linewidth(1.5)
ax[2].tick_params(which='minor', tickdir='in', length=3, top=True, right=True)
ax[2].tick_params(which='major', length=6)
ax[2].minorticks_on()
ax[2].xaxis.set_label_position('top')
ax[2].xaxis.label.set_color(axis_colors[1])
ax[2].set_xlabel(r'$Orbit\ type\ 9$', fontsize=35)
ax[2].set_ylim(-0.7, 2.4)
ax[2].axvline(450e-3, color='b', linestyle='dashed', linewidth=2.0, alpha=0.6)

plt.subplots_adjust(wspace=0, hspace=0)

plt.show()

os.chdir('/home/elismar/Documentos/Fisica/IC/GalMer/inflow')
