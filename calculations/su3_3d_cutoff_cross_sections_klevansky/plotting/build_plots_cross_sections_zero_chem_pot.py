import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import *
from common_utils.cross_section_data import *


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
data_folder = "su3_3d_cutoff_cross_sections_klevansky/data/"
plots_folder = "su3_3d_cutoff_cross_sections_klevansky/plots/"


####################################################################################################
# Cross section for uu->uu , ud->ud
print("Building plot: cross section for uu->uu , ud->ud (T[GeV]=0.215, zero chemical potential)")

# Load data from the file
filename = "crossSectionUUUU_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uuuu_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udud_T0215 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uuuu = 'black'
linestyle_uuuu = '-'
color_udud = 'red'
linestyle_udud = '-'

# Plot cross sections
ax.plot(data_uuuu_T0215.get_sqrt_center_of_mass_energy(), data_uuuu_T0215.get_cross_section(), label=r'$uu \rightarrow uu$', color=color_uuuu, linewidth=2, linestyle=linestyle_uuuu)
ax.plot(data_udud_T0215.get_sqrt_center_of_mass_energy(), data_udud_T0215.get_cross_section(), label=r'$ud \rightarrow ud$', color=color_udud, linewidth=2, linestyle=linestyle_udud)

# Plot minimum center of mass energy point
ax.plot(data_uuuu_T0215.get_sqrt_center_of_mass_energy()[0], data_uuuu_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_uuuu)
ax.plot(data_udud_T0215.get_sqrt_center_of_mass_energy()[0], data_udud_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_udud)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.215', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uuuu_udud_T0215_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uuuu_T0215, data_udud_T0215
del color_uuuu, color_udud
del linestyle_uuuu, linestyle_udud
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for uu->uu , ud->ud
print("Building plot: cross section for uu->uu , ud->ud (T[GeV]=0.250, zero chemical potential)")

# Load data from the file
filename = "crossSectionUUUU_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uuuu_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udud_T0250 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uuuu = 'black'
linestyle_uuuu = '-'
color_udud = 'red'
linestyle_udud = '-'

# Plot cross sections
ax.plot(data_uuuu_T0250.get_sqrt_center_of_mass_energy(), data_uuuu_T0250.get_cross_section(), label=r'$uu \rightarrow uu$', color=color_uuuu, linewidth=2, linestyle=linestyle_uuuu)
ax.plot(data_udud_T0250.get_sqrt_center_of_mass_energy(), data_udud_T0250.get_cross_section(), label=r'$ud \rightarrow ud$', color=color_udud, linewidth=2, linestyle=linestyle_udud)

# Plot minimum center of mass energy point
ax.plot(data_uuuu_T0250.get_sqrt_center_of_mass_energy()[0], data_uuuu_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_uuuu)
ax.plot(data_udud_T0250.get_sqrt_center_of_mass_energy()[0], data_udud_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_udud)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uuuu_udud_T0250_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uuuu_T0250, data_udud_T0250
del color_uuuu, color_udud
del linestyle_uuuu, linestyle_udud
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for us->us , ss->ss
print("Building plot: cross section for us->us , ss->ss (T[GeV]=0.215, zero chemical potential)")

# Load data from the file
filename = "crossSectionUSUS_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usus_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssss_T0215 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usus = 'black'
linestyle_usus = '-'
color_ssss = 'red'
linestyle_ssss = '-'

# Plot cross sections
ax.plot(data_usus_T0215.get_sqrt_center_of_mass_energy(), data_usus_T0215.get_cross_section(), label=r'$us \rightarrow us$', color=color_usus, linewidth=2, linestyle=linestyle_usus)
ax.plot(data_ssss_T0215.get_sqrt_center_of_mass_energy(), data_ssss_T0215.get_cross_section(), label=r'$ss \rightarrow ss$', color=color_ssss, linewidth=2, linestyle=linestyle_ssss)

# Plot minimum center of mass energy point
ax.plot(data_usus_T0215.get_sqrt_center_of_mass_energy()[0], data_usus_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_usus)
ax.plot(data_ssss_T0215.get_sqrt_center_of_mass_energy()[0], data_ssss_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_ssss)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.215', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usus_ssss_T0215_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usus_T0215, data_ssss_T0215
del color_usus, color_ssss
del linestyle_usus, linestyle_ssss
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for us->us , ss->ss
print("Building plot: cross section for us->us , ss->ss (T[GeV]=0.250, zero chemical potential)")

# Load data from the file
filename = "crossSectionUSUS_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usus_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssss_T0250 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usus = 'black'
linestyle_usus = '-'
color_ssss = 'red'
linestyle_ssss = '-'

# Plot cross sections
ax.plot(data_usus_T0250.get_sqrt_center_of_mass_energy(), data_usus_T0250.get_cross_section(), label=r'$us \rightarrow us$', color=color_usus, linewidth=2, linestyle=linestyle_usus)
ax.plot(data_ssss_T0250.get_sqrt_center_of_mass_energy(), data_ssss_T0250.get_cross_section(), label=r'$ss \rightarrow ss$', color=color_ssss, linewidth=2, linestyle=linestyle_ssss)

# Plot minimum center of mass energy point
ax.plot(data_usus_T0250.get_sqrt_center_of_mass_energy()[0], data_usus_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_usus)
ax.plot(data_ssss_T0250.get_sqrt_center_of_mass_energy()[0], data_ssss_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_ssss)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 6.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usus_ssss_T0250_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usus_T0250, data_ssss_T0250
del color_usus, color_ssss
del linestyle_usus, linestyle_ssss
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for udbar->udbar , uubar->uubar , uubar->ddbar , uubar->ssbar
print("Building plot: cross section for udbar->udbar , uubar->uubar, uubar->ddbar, uubar->ssbar (T[GeV]=0.215, zero chemical potential)")

# Load data from the file
filename = "crossSectionUDBARUDBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udbarudbar_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARUUBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubaruubar_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARDDBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarddbar_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARSSBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarssbar_T0215 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_udbarudbar = 'black'
linestyle_udbarudbar = '-'
color_uubaruubar = 'red'
linestyle_uubaruubar = '-'
color_uubarddbar = 'blue'
linestyle_uubarddbar = '-'
color_uubarssbar = 'lime'
linestyle_uubarssbar = '-'


# Plot cross sections
#ax.plot([], [], ' ', label=r'T[GeV]=0.215')
ax.plot(data_udbarudbar_T0215.get_sqrt_center_of_mass_energy(), data_udbarudbar_T0215.get_cross_section(), label=r'$u\bar{d} \rightarrow u\bar{d}$', color=color_udbarudbar, linewidth=2, linestyle=linestyle_udbarudbar)
ax.plot(data_uubaruubar_T0215.get_sqrt_center_of_mass_energy(), data_uubaruubar_T0215.get_cross_section(), label=r'$u\bar{u} \rightarrow u\bar{u}$', color=color_uubaruubar, linewidth=2, linestyle=linestyle_uubaruubar)
ax.plot(data_uubarddbar_T0215.get_sqrt_center_of_mass_energy(), data_uubarddbar_T0215.get_cross_section(), label=r'$u\bar{u} \rightarrow d\bar{d}$', color=color_uubarddbar, linewidth=2, linestyle=linestyle_uubarddbar)
ax.plot(data_uubarssbar_T0215.get_sqrt_center_of_mass_energy(), data_uubarssbar_T0215.get_cross_section(), label=r'$u\bar{u} \rightarrow s\bar{s}$', color=color_uubarssbar, linewidth=2, linestyle=linestyle_uubarssbar)

# Plot minimum center of mass energy point
ax.plot(data_udbarudbar_T0215.get_sqrt_center_of_mass_energy()[0], data_udbarudbar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_udbarudbar)
ax.plot(data_uubaruubar_T0215.get_sqrt_center_of_mass_energy()[0], data_uubaruubar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_uubaruubar)
ax.plot(data_uubaruubar_T0215.get_sqrt_center_of_mass_energy()[0], data_uubarddbar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_uubarddbar)
ax.plot(data_uubarssbar_T0215.get_sqrt_center_of_mass_energy()[0], data_uubarssbar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_uubarssbar)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=14, frameon=False, title='T[GeV]=0.215', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 50.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_udbarudbar_uubaruubar_uubarddbar_uubarssbar_T0215_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_udbarudbar_T0215, data_uubaruubar_T0215, data_uubarddbar_T0215, data_uubarssbar_T0215
del color_udbarudbar, color_uubaruubar, color_uubarddbar, color_uubarssbar
del linestyle_udbarudbar, linestyle_uubaruubar, linestyle_uubarddbar, linestyle_uubarssbar
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for udbar->udbar , uubar->uubar , uubar->ddbar , uubar->ssbar
print("Building plot: cross section for udbar->udbar , uubar->uubar, uubar->ddbar, uubar->ssbar (T[GeV]=0.250, zero chemical potential)")

# Load data from the file
filename = "crossSectionUDBARUDBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udbarudbar_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARUUBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubaruubar_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARDDBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarddbar_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBARSSBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarssbar_T0250 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_udbarudbar = 'black'
linestyle_udbarudbar = '-'
color_uubaruubar = 'red'
linestyle_uubaruubar = '-'
color_uubarddbar = 'blue'
linestyle_uubarddbar = '-'
color_uubarssbar = 'lime'
linestyle_uubarssbar = '-'


# Plot cross sections
ax.plot(data_udbarudbar_T0250.get_sqrt_center_of_mass_energy(), data_udbarudbar_T0250.get_cross_section(), label=r'$u\bar{d} \rightarrow u\bar{d}$', color=color_udbarudbar, linewidth=2, linestyle=linestyle_udbarudbar)
ax.plot(data_uubaruubar_T0250.get_sqrt_center_of_mass_energy(), data_uubaruubar_T0250.get_cross_section(), label=r'$u\bar{u} \rightarrow u\bar{u}$', color=color_uubaruubar, linewidth=2, linestyle=linestyle_uubaruubar)
ax.plot(data_uubarddbar_T0250.get_sqrt_center_of_mass_energy(), data_uubarddbar_T0250.get_cross_section(), label=r'$u\bar{u} \rightarrow d\bar{d}$', color=color_uubarddbar, linewidth=2, linestyle=linestyle_uubarddbar)
ax.plot(data_uubarssbar_T0250.get_sqrt_center_of_mass_energy(), data_uubarssbar_T0250.get_cross_section(), label=r'$u\bar{u} \rightarrow s\bar{s}$', color=color_uubarssbar, linewidth=2, linestyle=linestyle_uubarssbar)

# Plot minimum center of mass energy point
ax.plot(data_udbarudbar_T0250.get_sqrt_center_of_mass_energy()[0], data_udbarudbar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_udbarudbar)
ax.plot(data_uubaruubar_T0250.get_sqrt_center_of_mass_energy()[0], data_uubaruubar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_uubaruubar)
ax.plot(data_uubaruubar_T0250.get_sqrt_center_of_mass_energy()[0], data_uubarddbar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_uubarddbar)
ax.plot(data_uubarssbar_T0250.get_sqrt_center_of_mass_energy()[0], data_uubarssbar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_uubarssbar)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=14, frameon=False, title='T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 50.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_udbarudbar_uubaruubar_uubarddbar_uubarssbar_T0250_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_udbarudbar_T0250, data_uubaruubar_T0250, data_uubarddbar_T0250, data_uubarssbar_T0250
del color_udbarudbar, color_uubaruubar, color_uubarddbar, color_uubarssbar
del linestyle_udbarudbar, linestyle_uubaruubar, linestyle_uubarddbar, linestyle_uubarssbar
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for usbar->usbar , ssbar->uubar , ssbar->ssbar 
print("Building plot: cross section for usbar->usbar , ssbar->uubar, ssbar->ssbar (T[GeV]=0.215, zero chemical potential)")

# Load data from the file
filename = "crossSectionUSBARUSBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usbarusbar_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBARUUBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbaruubar_T0215 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBARSSBAR_T0.215000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbarssbar_T0215 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usbarusbar = 'black'
linestyle_usbarusbar = '-'
color_ssbaruubar = 'red'
linestyle_ssbaruubar = '-'
color_ssbarssbar = 'blue'
linestyle_ssbarssbar = '-'


# Plot cross sections
ax.plot(data_usbarusbar_T0215.get_sqrt_center_of_mass_energy(), data_usbarusbar_T0215.get_cross_section(), label=r'$u\bar{s} \rightarrow u\bar{s}$', color=color_usbarusbar, linewidth=2, linestyle=linestyle_usbarusbar)
ax.plot(data_ssbaruubar_T0215.get_sqrt_center_of_mass_energy(), data_ssbaruubar_T0215.get_cross_section(), label=r'$s\bar{s} \rightarrow u\bar{u}$', color=color_ssbaruubar, linewidth=2, linestyle=linestyle_ssbaruubar)
ax.plot(data_ssbarssbar_T0215.get_sqrt_center_of_mass_energy(), data_ssbarssbar_T0215.get_cross_section(), label=r'$s\bar{s} \rightarrow s\bar{s}$', color=color_ssbarssbar, linewidth=2, linestyle=linestyle_ssbarssbar)

# Plot minimum center of mass energy point
ax.plot(data_usbarusbar_T0215.get_sqrt_center_of_mass_energy()[0], data_usbarusbar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_usbarusbar)
ax.plot(data_ssbaruubar_T0215.get_sqrt_center_of_mass_energy()[0], data_ssbaruubar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_ssbaruubar)
ax.plot(data_ssbaruubar_T0215.get_sqrt_center_of_mass_energy()[0], data_ssbarssbar_T0215.get_cross_section()[0], marker='o', markersize=6, color=color_ssbarssbar)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.215', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 25.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usbarusbar_ssbaruubar_ssbarssbar_T0215_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usbarusbar_T0215, data_ssbaruubar_T0215, data_ssbarssbar_T0215
del color_usbarusbar, color_ssbaruubar, color_ssbarssbar
del linestyle_usbarusbar, linestyle_ssbaruubar, linestyle_ssbarssbar
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for usbar->usbar , ssbar->uubar , ssbar->ssbar 
print("Building plot: cross section for usbar->usbar , ssbar->uubar, ssbar->ssbar (T[GeV]=0.250, zero chemical potential)")

# Load data from the file
filename = "crossSectionUSBARUSBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usbarusbar_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBARUUBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbaruubar_T0250 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBARSSBAR_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbarssbar_T0250 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usbarusbar = 'black'
linestyle_usbarusbar = '-'
color_ssbaruubar = 'red'
linestyle_ssbaruubar = '-'
color_ssbarssbar = 'blue'
linestyle_ssbarssbar = '-'


# Plot cross sections
ax.plot(data_usbarusbar_T0250.get_sqrt_center_of_mass_energy(), data_usbarusbar_T0250.get_cross_section(), label=r'$u\bar{s} \rightarrow u\bar{s}$', color=color_usbarusbar, linewidth=2, linestyle=linestyle_usbarusbar)
ax.plot(data_ssbaruubar_T0250.get_sqrt_center_of_mass_energy(), data_ssbaruubar_T0250.get_cross_section(), label=r'$s\bar{s} \rightarrow u\bar{u}$', color=color_ssbaruubar, linewidth=2, linestyle=linestyle_ssbaruubar)
ax.plot(data_ssbarssbar_T0250.get_sqrt_center_of_mass_energy(), data_ssbarssbar_T0250.get_cross_section(), label=r'$s\bar{s} \rightarrow s\bar{s}$', color=color_ssbarssbar, linewidth=2, linestyle=linestyle_ssbarssbar)

# Plot minimum center of mass energy point
ax.plot(data_usbarusbar_T0250.get_sqrt_center_of_mass_energy()[0], data_usbarusbar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_usbarusbar)
ax.plot(data_ssbaruubar_T0250.get_sqrt_center_of_mass_energy()[0], data_ssbaruubar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_ssbaruubar)
ax.plot(data_ssbaruubar_T0250.get_sqrt_center_of_mass_energy()[0], data_ssbarssbar_T0250.get_cross_section()[0], marker='o', markersize=6, color=color_ssbarssbar)

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title='T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 25.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usbarusbar_ssbaruubar_ssbarssbar_T0250_CP0.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usbarusbar_T0250, data_ssbaruubar_T0250, data_ssbarssbar_T0250
del color_usbarusbar, color_ssbaruubar, color_ssbarssbar
del linestyle_usbarusbar, linestyle_ssbaruubar, linestyle_ssbarssbar
del fig, ax, xmin, xmax, ymin, ymax
