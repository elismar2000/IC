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

cmap = plt.get_cmap('twilight')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig, axs = plt.subplots(4, 1, figsize=(8, 6), sharey=False, sharex=True)

time = np.arange(1, 27, 1) * 0.05 * 0.98 #Time in Gyr

#============================
#Radius = 1.0 Kpc
#============================

gm1 = np.array([ j['gal1_r10'] for j in inflow_table.iterrows()])
gm2 = np.array([ j['gal2_r10'] for j in inflow_table.iterrows()])
gm = np.array([ j['gal_r10'] for j in inflow_table.iterrows()])

n1 = np.array([ j['gal1_r10'] for j in part_num_table.iterrows()])
n2 = np.array([ j['gal2_r10'] for j in part_num_table.iterrows()])
n = np.array([ j['gal_r10'] for j in part_num_table.iterrows()])

gm1 = gm1[0:27]
gm2 = gm2[0:27]
gm = gm[0:27]

n1 = n1[0:27]
n2 = n2[0:27]
n = n[0:27]

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

err1 = (np.sqrt(n1) * gm1 / n1) / age_factor
err2 = (np.sqrt(n2) * gm2 / n2) / age_factor
err = (np.sqrt(n) * gm / n) / age_factor

flux1 = np.diff(gm1) / age_factor
flux2 = np.diff(gm2) / age_factor
flux = np.diff(gm) / age_factor

axs[0].plot(time, (flux1 * mass_factor) / 10, '-', color=colors[3], linewidth=1.0)
axs[0].plot(time, (flux2 * mass_factor) / 10, '-', color=colors[1], linewidth=1.0)
axs[0].plot(time, (flux * mass_factor) / 10, '-', color=colors[5], linewidth=1.0)

#10 dividindo ali é pra tirar a notação científica do eixo y
axs[0].fill_between(time, (mass_factor/10) * (flux1 - err1[1:]), (mass_factor/10) * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
axs[0].fill_between(time, (mass_factor/10) * (flux2 - err2[1:]), (mass_factor/10) * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
axs[0].fill_between(time, (mass_factor/10) * (flux - err[1:]), (mass_factor/10) * (flux + err[1:]), alpha=0.3, color=colors[5])

axs[0].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[0].tick_params(which='minor', length=1)
axs[0].tick_params(which='major', length=2)
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[0].minorticks_on()

axs[0].set_ylabel(r'$Flux\ [10 \times \frac{M\odot}{yr}]$', fontsize=10, labelpad=15)

axs[0].text(0.0, -0.7, 'r = 1.0 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.5 Kpc
#============================

gm1 = np.array([ j['gal1_r05'] for j in inflow_table.iterrows()])
gm2 = np.array([ j['gal2_r05'] for j in inflow_table.iterrows()])
gm = np.array([ j['gal_r05'] for j in inflow_table.iterrows()])

n1 = np.array([ j['gal1_r05'] for j in part_num_table.iterrows()])
n2 = np.array([ j['gal2_r05'] for j in part_num_table.iterrows()])
n = np.array([ j['gal_r05'] for j in part_num_table.iterrows()])

gm1 = gm1[0:27]
gm2 = gm2[0:27]
gm = gm[0:27]

n1 = n1[0:27]
n2 = n2[0:27]
n = n[0:27]

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

err1 = (np.sqrt(n1) * gm1 / n1) / age_factor
err2 = (np.sqrt(n2) * gm2 / n2) / age_factor
err = (np.sqrt(n) * gm / n) / age_factor

flux1 = np.diff(gm1) / age_factor
flux2 = np.diff(gm2) / age_factor
flux = np.diff(gm) / age_factor

axs[1].plot(time, (flux1 * mass_factor) / 10, '-', color=colors[3], linewidth=1.0)
axs[1].plot(time, (flux2 * mass_factor) / 10, '-', color=colors[1], linewidth=1.0)
axs[1].plot(time, (flux * mass_factor) / 10, '-', color=colors[5], linewidth=1.0)

#10 dividindo ali é pra tirar a notação científica do eixo y
axs[1].fill_between(time, (mass_factor/10) * (flux1 - err1[1:]), (mass_factor/10) * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
axs[1].fill_between(time, (mass_factor/10) * (flux2 - err2[1:]), (mass_factor/10) * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
axs[1].fill_between(time, (mass_factor/10) * (flux - err[1:]), (mass_factor/10) * (flux + err[1:]), alpha=0.3, color=colors[5])

axs[1].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[1].tick_params(which='minor', length=1)
axs[1].tick_params(which='major', length=2)
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[1].minorticks_on()

axs[1].set_ylabel(r'$Flux\ [10 \times \frac{M\odot}{yr}]$', fontsize=10)

axs[1].text(0.0, -0.8, 'r = 0.5 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.1 Kpc
#============================

gm1 = np.array([ j['gal1_r01'] for j in inflow_table.iterrows()])
gm2 = np.array([ j['gal2_r01'] for j in inflow_table.iterrows()])
gm = np.array([ j['gal_r01'] for j in inflow_table.iterrows()])

n1 = np.array([ j['gal1_r01'] for j in part_num_table.iterrows()])
n2 = np.array([ j['gal2_r01'] for j in part_num_table.iterrows()])
n = np.array([ j['gal_r01'] for j in part_num_table.iterrows()])

gm1 = gm1[0:27]
gm2 = gm2[0:27]
gm = gm[0:27]

n1 = n1[0:27]
n2 = n2[0:27]
n = n[0:27]

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

err1 = (np.sqrt(n1) * gm1 / n1) / age_factor
err2 = (np.sqrt(n2) * gm2 / n2) / age_factor
err = (np.sqrt(n) * gm / n) / age_factor

flux1 = np.diff(gm1) / age_factor
flux2 = np.diff(gm2) / age_factor
flux = np.diff(gm) / age_factor

axs[2].plot(time, (flux1 * mass_factor) / 10, '-', color=colors[3], linewidth=1.5)
axs[2].plot(time, (flux2 * mass_factor) / 10, '-', color=colors[1], linewidth=1.5)
axs[2].plot(time, (flux * mass_factor) / 10, '-', color=colors[5], linewidth=1.5)

#10 dividindo ali é pra tirar a notação científica do eixo y
axs[2].fill_between(time, (mass_factor/10) * (flux1 - err1[1:]), (mass_factor/10) * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
axs[2].fill_between(time, (mass_factor/10) * (flux2 - err2[1:]), (mass_factor/10) * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
axs[2].fill_between(time, (mass_factor/10) * (flux - err[1:]), (mass_factor/10) * (flux + err[1:]), alpha=0.3, color=colors[5])

axs[2].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[2].tick_params(which='minor', length=1)
axs[2].tick_params(which='major', length=2)
axs[2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[2].minorticks_on()

axs[2].set_ylabel(r'$Flux\ [10 \times \frac{M\odot}{yr}]$', fontsize=10, labelpad=15)

axs[2].text(0.0, -0.45, 'r = 0.1 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#Radius = 0.01 Kpc
#============================

gm1 = np.array([ j['gal1_r001'] for j in inflow_table.iterrows()])
gm2 = np.array([ j['gal2_r001'] for j in inflow_table.iterrows()])
gm = np.array([ j['gal_r001'] for j in inflow_table.iterrows()])

n1 = np.array([ j['gal1_r001'] for j in part_num_table.iterrows()])
n2 = np.array([ j['gal2_r001'] for j in part_num_table.iterrows()])
n = np.array([ j['gal_r001'] for j in part_num_table.iterrows()])

gm1 = gm1[0:27]
gm2 = gm2[0:27]
gm = gm[0:27]

n1 = n1[0:27]
n2 = n2[0:27]
n = n[0:27]

age_factor = 50e+6 * 0.98 #Time between each snapshot of the simulation in years
mass_factor = 1e+10

err1 = (np.sqrt(n1) * gm1 / n1) / age_factor
err2 = (np.sqrt(n2) * gm2 / n2) / age_factor
err = (np.sqrt(n) * gm / n) / age_factor

flux1 = np.diff(gm1) / age_factor
flux2 = np.diff(gm2) / age_factor
flux = np.diff(gm) / age_factor

axs[3].plot(time, (flux1 * mass_factor) / 1e-1, '-', color=colors[3], linewidth=1.5)
axs[3].plot(time, (flux2 * mass_factor) / 1e-1, '-', color=colors[1], linewidth=1.5)
axs[3].plot(time, (flux * mass_factor) / 1e-1, '-', color=colors[5], linewidth=1.5)

#10 dividindo ali é pra tirar a notação científica do eixo y
axs[3].fill_between(time, (mass_factor/1e-1) * (flux1 - err1[1:]), (mass_factor/1e-1) * (flux1 + err1[1:]), alpha=0.3, color=colors[3])
axs[3].fill_between(time, (mass_factor/1e-1) * (flux2 - err2[1:]), (mass_factor/1e-1) * (flux2 + err2[1:]), alpha=0.3, color=colors[1])
axs[3].fill_between(time, (mass_factor/1e-1) * (flux - err[1:]), (mass_factor/1e-1) * (flux + err[1:]), alpha=0.3, color=colors[5])

axs[3].tick_params(axis='both', which='both', tickdir='in', labelsize=15,
    top=True, right=True, width=1.5)
axs[3].tick_params(which='minor', length=1)
axs[3].tick_params(which='major', length=2)
axs[3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axs[3].minorticks_on()

axs[3].set_ylabel(r'$Flux\ [10^{-1} \times \frac{M\odot}{yr}]$', fontsize=10, labelpad=15)
axs[3].set_xlabel(r'$Time\ [Gyr]$', fontsize=20)

axs[3].text(0.0, -0.35, 'r = 0.01 Kpc', fontsize=10, bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)))

#============================
#For GalMer
#============================

f = 2.25e+9

time1 = np.arange(1, 17, 1) * 0.05 #Gyr
time2 = np.arange(17, 27, 1) * 0.05 #Gyr

#============================
#Radius 1.0 Kpc
#============================
gm1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm1.orbit1.min_radius1.0.npy')
gm2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm2.orbit1.min_radius1.0.npy')
gm = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm.orbit1.min_radius1.0.npy')

n1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n1.orbit1.min_radius1.0.npy')
n2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n2.orbit1.min_radius1.0.npy')
n = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n.orbit1.min_radius1.0.npy')

gm1 = gm1[0:17]
gm2 = gm2[0:17]
gm = gm[0:10]

n1 = n1[0:17]
n2 = n2[0:17]
n = n[0:10]

flux1 = np.diff(gm1) / 50e+6
flux2 = np.diff(gm2) / 50e+6
flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
err = (np.sqrt(n) * gm / n) / 50e+6

axs[0].plot(time1, (flux1 * f) / 10, '--', color=colors[3], linewidth=1.0, alpha=0.4)
axs[0].plot(time1, (flux2 * f) / 10, '--', color=colors[1], linewidth=1.0, alpha=0.4)
axs[0].plot(time2, (flux * f) / 10, '--', color=colors[5], linewidth=1.0, alpha=0.4)

axs[0].fill_between(time1, (f / 10) * (flux1 - err1[1:]), (f / 10) * (flux1 + err1[1:]), alpha=0.1, color=colors[3])
axs[0].fill_between(time1, (f / 10) * (flux2 - err2[1:]), (f / 10) * (flux2 + err2[1:]), alpha=0.1, color=colors[1])
axs[0].fill_between(time2, (f / 10) * (flux - err), (f / 10) * (flux + err), alpha=0.1, color=colors[5])

#============================
#Radius 0.5 Kpc
#============================

gm1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm1.orbit1.min_radius0.5.npy')
gm2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm2.orbit1.min_radius0.5.npy')
gm = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm.orbit1.min_radius0.5.npy')

n1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n1.orbit1.min_radius0.5.npy')
n2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n2.orbit1.min_radius0.5.npy')
n = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n.orbit1.min_radius0.5.npy')

gm1 = gm1[0:17]
gm2 = gm2[0:17]
gm = gm[0:10]

n1 = n1[0:17]
n2 = n2[0:17]
n = n[0:10]

flux1 = np.diff(gm1) / 50e+6
flux2 = np.diff(gm2) / 50e+6
flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
err = (np.sqrt(n) * gm / n) / 50e+6

axs[1].plot(time1, (flux1 * f)/10, '--', color=colors[3], linewidth=1.5, alpha=0.4)
axs[1].plot(time1, (flux2 * f)/10, '--', color=colors[1], linewidth=1.5, alpha=0.4)
axs[1].plot(time2, (flux * f)/10, '--', color=colors[5], linewidth=1.5, alpha=0.4)

axs[1].fill_between(time1, f * (flux1 - err1[1:])/10, f * (flux1 + err1[1:])/10, alpha=0.1, color=colors[3])
axs[1].fill_between(time1, f * (flux2 - err2[1:])/10, f * (flux2 + err2[1:])/10, alpha=0.1, color=colors[1])
axs[1].fill_between(time2, f * (flux - err)/10, f * (flux + err)/10, alpha=0.1, color=colors[5])

#============================
#Radius 0.1 Kpc
#============================
gm1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm1.orbit1.min_radius0.1.npy')
gm2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm2.orbit1.min_radius0.1.npy')
gm = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm.orbit1.min_radius0.1.npy')

n1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n1.orbit1.min_radius0.1.npy')
n2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n2.orbit1.min_radius0.1.npy')
n = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n.orbit1.min_radius0.1.npy')

gm1 = gm1[0:17]
gm2 = gm2[0:17]
gm = gm[0:10]

n1 = n1[0:17]
n2 = n2[0:17]
n = n[0:10]

flux1 = np.diff(gm1) / 50e+6
flux2 = np.diff(gm2) / 50e+6
flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
err = (np.sqrt(n) * gm / n) / 50e+6

axs[2].plot(time1, (flux1 * f)/10, '--', color=colors[3], linewidth=1.5, alpha=0.4)
axs[2].plot(time1, (flux2 * f)/10, '--', color=colors[1], linewidth=1.5, alpha=0.4)
axs[2].plot(time2, (flux * f)/10, '--', color=colors[5], linewidth=1.5, alpha=0.4)

axs[2].fill_between(time1, f * (flux1 - err1[1:])/10, f * (flux1 + err1[1:])/10, alpha=0.1, color=colors[3])
axs[2].fill_between(time1, f * (flux2 - err2[1:])/10, f * (flux2 + err2[1:])/10, alpha=0.1, color=colors[1])
axs[2].fill_between(time2, f * (flux - err)/10, f * (flux + err)/10, alpha=0.1, color=colors[5])

#============================
#Radius 0.01 Kpc
#============================

gm1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm1.orbit1.min_radius0.01.npy')
gm2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm2.orbit1.min_radius0.01.npy')
gm = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/inflow_outputs_galmer/m_gm.orbit1.min_radius0.01.npy')

n1 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n1.orbit1.min_radius0.01.npy')
n2 = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n2.orbit1.min_radius0.01.npy')
n = np.load('/home/elismar/Documentos/Fisica/IC/GalMer/inflow/part_num_outputs_galmer/n.orbit1.min_radius0.01.npy')

gm1 = gm1[0:17]
gm2 = gm2[0:17]
gm = gm[0:10]

n1 = n1[0:17]
n2 = n2[0:17]
n = n[0:10]

flux1 = np.diff(gm1) / 50e+6
flux2 = np.diff(gm2) / 50e+6
flux = np.concatenate(([gm[0] - (gm1[-1] + gm2[-1])], np.diff(gm))) / 50e+6

err1 = (np.sqrt(n1) * gm1 / n1) / 50e+6
err2 = (np.sqrt(n2) * gm2 / n2) / 50e+6
err = (np.sqrt(n) * gm / n) / 50e+6

axs[3].plot(time1, (flux1 * f)/1e-1, '--', color=colors[3], linewidth=1.5, alpha=0.4)
axs[3].plot(time1, (flux2 * f)/1e-1, '--', color=colors[1], linewidth=1.5, alpha=0.4)
axs[3].plot(time2, (flux * f)/1e-1, '--', color=colors[5], linewidth=1.5, alpha=0.4)

axs[3].fill_between(time1, f * (flux1 - err1[1:])/1e-1, f * (flux1 + err1[1:])/1e-1, alpha=0.1, color=colors[3])
axs[3].fill_between(time1, f * (flux2 - err2[1:])/1e-1, f * (flux2 + err2[1:])/1e-1, alpha=0.1, color=colors[1])
axs[3].fill_between(time2, f * (flux - err)/1e-1, f * (flux + err)/1e-1, alpha=0.1, color=colors[5])

#============================

axs[3].plot([0, 0], [0, 0], '--', color='k', label='GalMer')
axs[3].plot([0, 0], [0, 0], '-', color='k', label='GADGET-3')
plt.legend(loc=0, fontsize=10, fancybox=True, framealpha=0.2)

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('inflow_gadget_galmer.png')
#plt.show()

#Lembrar de fazer um loop da próxima vez que tiver que plotar algo repetitivo assim!!!!!
