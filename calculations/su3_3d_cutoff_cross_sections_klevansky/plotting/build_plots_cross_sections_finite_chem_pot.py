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
# Cross section for uu->uu
print("Building plot: cross section for uu->uu (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUUUU_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uuuu_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUUU_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_uuuu_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUUU_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_uuuu_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUUU_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_uuuu_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUUU_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_uuuu_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUUU_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_uuuu_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uuuu_T0250_Cp0000 = 'black'
linestyle_uuuu_T0250_Cp0000 = '-'
color_uuuu_T0250_Cp0100 = 'red'
linestyle_uuuu_T0250_Cp0100 = '-'
color_uuuu_T0250_Cp0200 = 'blue'
linestyle_uuuu_T0250_Cp0200 = '-'
color_uuuu_T0250_Cp0300 = 'lime'
linestyle_uuuu_T0250_Cp0300 = '-'
color_uuuu_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_uuuu_T0250_Cp0400 = '-'
color_uuuu_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_uuuu_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_uuuu_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_uuuu_T0250_Cp0000, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0000
    )
ax.plot(
    data_uuuu_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_uuuu_T0250_Cp0100, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0100
    )
ax.plot(
    data_uuuu_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_uuuu_T0250_Cp0200, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0200
    )
ax.plot(
    data_uuuu_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_uuuu_T0250_Cp0300, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0300
    )
ax.plot(
    data_uuuu_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_uuuu_T0250_Cp0400, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0400
    )
ax.plot(
    data_uuuu_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_uuuu_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_uuuu_T0250_Cp0500, linewidth=2, linestyle=linestyle_uuuu_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_uuuu_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0000
    )
ax.plot(
    data_uuuu_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0100
    )
ax.plot(
    data_uuuu_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0200
    )
ax.plot(
    data_uuuu_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0300
    )
ax.plot(
    data_uuuu_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0400
    )
ax.plot(
    data_uuuu_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_uuuu_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uuuu_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{uu \rightarrow uu} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
#plt.legend(loc='upper left', fontsize=14, frameon=False)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

# auxH = 0.065; auxX = 0.05; auxY = 0.25
# texts = [
#     r'$uu \rightarrow uu$',
#     r'T[GeV]=0.250'
# ]
# add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uuuu_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uuuu_T0250_Cp0000, data_uuuu_T0250_Cp0100, data_uuuu_T0250_Cp0200, data_uuuu_T0250_Cp0300, data_uuuu_T0250_Cp0400, data_uuuu_T0250_Cp0500
del color_uuuu_T0250_Cp0000, color_uuuu_T0250_Cp0100, color_uuuu_T0250_Cp0200, color_uuuu_T0250_Cp0300, color_uuuu_T0250_Cp0400, color_uuuu_T0250_Cp0500
del linestyle_uuuu_T0250_Cp0000, linestyle_uuuu_T0250_Cp0100, linestyle_uuuu_T0250_Cp0200, linestyle_uuuu_T0250_Cp0300, linestyle_uuuu_T0250_Cp0400, linestyle_uuuu_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ubarubar->ubarubar
print("Building plot: cross section for ubarubar->ubarubar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ubarubarubarubar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ubarubarubarubar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ubarubarubarubar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ubarubarubarubar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ubarubarubarubar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarUBarUBarUBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ubarubarubarubar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ubarubarubarubar_T0250_Cp0000 = 'black'
linestyle_ubarubarubarubar_T0250_Cp0000 = '-'
color_ubarubarubarubar_T0250_Cp0100 = 'red'
linestyle_ubarubarubarubar_T0250_Cp0100 = '-'
color_ubarubarubarubar_T0250_Cp0200 = 'blue'
linestyle_ubarubarubarubar_T0250_Cp0200 = '-'
color_ubarubarubarubar_T0250_Cp0300 = 'lime'
linestyle_ubarubarubarubar_T0250_Cp0300 = '-'
color_ubarubarubarubar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ubarubarubarubar_T0250_Cp0400 = '-'
color_ubarubarubarubar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ubarubarubarubar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ubarubarubarubar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ubarubarubarubar_T0250_Cp0000, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0000
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ubarubarubarubar_T0250_Cp0100, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0100
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ubarubarubarubar_T0250_Cp0200, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0200
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ubarubarubarubar_T0250_Cp0300, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0300
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ubarubarubarubar_T0250_Cp0400, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0400
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ubarubarubarubar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ubarubarubarubar_T0250_Cp0500, linewidth=2, linestyle=linestyle_ubarubarubarubar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ubarubarubarubar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0000
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0100
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0200
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0300
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0400
    )
ax.plot(
    data_ubarubarubarubar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ubarubarubarubar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarubarubarubar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\bar{u}\bar{u} \rightarrow \bar{u}\bar{u}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ubarubarubarubar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ubarubarubarubar_T0250_Cp0000, data_ubarubarubarubar_T0250_Cp0100, data_ubarubarubarubar_T0250_Cp0200, data_ubarubarubarubar_T0250_Cp0300, data_ubarubarubarubar_T0250_Cp0400, data_ubarubarubarubar_T0250_Cp0500
del color_ubarubarubarubar_T0250_Cp0000, color_ubarubarubarubar_T0250_Cp0100, color_ubarubarubarubar_T0250_Cp0200, color_ubarubarubarubar_T0250_Cp0300, color_ubarubarubarubar_T0250_Cp0400, color_ubarubarubarubar_T0250_Cp0500
del linestyle_ubarubarubarubar_T0250_Cp0000, linestyle_ubarubarubarubar_T0250_Cp0100, linestyle_ubarubarubarubar_T0250_Cp0200, linestyle_ubarubarubarubar_T0250_Cp0300, linestyle_ubarubarubarubar_T0250_Cp0400, linestyle_ubarubarubarubar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ud->ud
print("Building plot: cross section for ud->ud (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUDUD_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udud_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_udud_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_udud_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_udud_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_udud_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDUD_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_udud_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_udud_T0250_Cp0000 = 'black'
linestyle_udud_T0250_Cp0000 = '-'
color_udud_T0250_Cp0100 = 'red'
linestyle_udud_T0250_Cp0100 = '-'
color_udud_T0250_Cp0200 = 'blue'
linestyle_udud_T0250_Cp0200 = '-'
color_udud_T0250_Cp0300 = 'lime'
linestyle_udud_T0250_Cp0300 = '-'
color_udud_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_udud_T0250_Cp0400 = '-'
color_udud_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_udud_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_udud_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_udud_T0250_Cp0000, linewidth=2, linestyle=linestyle_udud_T0250_Cp0000
    )
ax.plot(
    data_udud_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_udud_T0250_Cp0100, linewidth=2, linestyle=linestyle_udud_T0250_Cp0100
    )
ax.plot(
    data_udud_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_udud_T0250_Cp0200, linewidth=2, linestyle=linestyle_udud_T0250_Cp0200
    )
ax.plot(
    data_udud_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_udud_T0250_Cp0300, linewidth=2, linestyle=linestyle_udud_T0250_Cp0300
    )
ax.plot(
    data_udud_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_udud_T0250_Cp0400, linewidth=2, linestyle=linestyle_udud_T0250_Cp0400
    )
ax.plot(
    data_udud_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_udud_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_udud_T0250_Cp0500, linewidth=2, linestyle=linestyle_udud_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_udud_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0000
    )
ax.plot(
    data_udud_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0100
    )
ax.plot(
    data_udud_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0200
    )
ax.plot(
    data_udud_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0300
    )
ax.plot(
    data_udud_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0400
    )
ax.plot(
    data_udud_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_udud_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udud_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{ud \rightarrow ud} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_udud_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_udud_T0250_Cp0000, data_udud_T0250_Cp0100, data_udud_T0250_Cp0200, data_udud_T0250_Cp0300, data_udud_T0250_Cp0400, data_udud_T0250_Cp0500
del color_udud_T0250_Cp0000, color_udud_T0250_Cp0100, color_udud_T0250_Cp0200, color_udud_T0250_Cp0300, color_udud_T0250_Cp0400, color_udud_T0250_Cp0500
del linestyle_udud_T0250_Cp0000, linestyle_udud_T0250_Cp0100, linestyle_udud_T0250_Cp0200, linestyle_udud_T0250_Cp0300, linestyle_udud_T0250_Cp0400, linestyle_udud_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ubardbar->ubardbar
print("Building plot: cross section for ubardbar->ubardbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ubardbarubardbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ubardbarubardbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ubardbarubardbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ubardbarubardbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ubardbarubardbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarDBarUBarDBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ubardbarubardbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ubardbarubardbar_T0250_Cp0000 = 'black'
linestyle_ubardbarubardbar_T0250_Cp0000 = '-'
color_ubardbarubardbar_T0250_Cp0100 = 'red'
linestyle_ubardbarubardbar_T0250_Cp0100 = '-'
color_ubardbarubardbar_T0250_Cp0200 = 'blue'
linestyle_ubardbarubardbar_T0250_Cp0200 = '-'
color_ubardbarubardbar_T0250_Cp0300 = 'lime'
linestyle_ubardbarubardbar_T0250_Cp0300 = '-'
color_ubardbarubardbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ubardbarubardbar_T0250_Cp0400 = '-'
color_ubardbarubardbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ubardbarubardbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ubardbarubardbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ubardbarubardbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0000
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ubardbarubardbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0100
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ubardbarubardbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0200
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ubardbarubardbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0300
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ubardbarubardbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0400
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ubardbarubardbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ubardbarubardbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_ubardbarubardbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ubardbarubardbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0000
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0100
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0200
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0300
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0400
    )
ax.plot(
    data_ubardbarubardbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ubardbarubardbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubardbarubardbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\bar{u}\bar{d} \rightarrow \bar{u}\bar{d}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ubardbarubardbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ubardbarubardbar_T0250_Cp0000, data_ubardbarubardbar_T0250_Cp0100, data_ubardbarubardbar_T0250_Cp0200, data_ubardbarubardbar_T0250_Cp0300, data_ubardbarubardbar_T0250_Cp0400, data_ubardbarubardbar_T0250_Cp0500
del color_ubardbarubardbar_T0250_Cp0000, color_ubardbarubardbar_T0250_Cp0100, color_ubardbarubardbar_T0250_Cp0200, color_ubardbarubardbar_T0250_Cp0300, color_ubardbarubardbar_T0250_Cp0400, color_ubardbarubardbar_T0250_Cp0500
del linestyle_ubardbarubardbar_T0250_Cp0000, linestyle_ubardbarubardbar_T0250_Cp0100, linestyle_ubardbarubardbar_T0250_Cp0200, linestyle_ubardbarubardbar_T0250_Cp0300, linestyle_ubardbarubardbar_T0250_Cp0400, linestyle_ubardbarubardbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for us->us
print("Building plot: cross section for us->us (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUSUS_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usus_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSUS_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_usus_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSUS_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_usus_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSUS_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_usus_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSUS_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_usus_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSUS_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_usus_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usus_T0250_Cp0000 = 'black'
linestyle_usus_T0250_Cp0000 = '-'
color_usus_T0250_Cp0100 = 'red'
linestyle_usus_T0250_Cp0100 = '-'
color_usus_T0250_Cp0200 = 'blue'
linestyle_usus_T0250_Cp0200 = '-'
color_usus_T0250_Cp0300 = 'lime'
linestyle_usus_T0250_Cp0300 = '-'
color_usus_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_usus_T0250_Cp0400 = '-'
color_usus_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_usus_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_usus_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_usus_T0250_Cp0000, linewidth=2, linestyle=linestyle_usus_T0250_Cp0000
    )
ax.plot(
    data_usus_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_usus_T0250_Cp0100, linewidth=2, linestyle=linestyle_usus_T0250_Cp0100
    )
ax.plot(
    data_usus_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_usus_T0250_Cp0200, linewidth=2, linestyle=linestyle_usus_T0250_Cp0200
    )
ax.plot(
    data_usus_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_usus_T0250_Cp0300, linewidth=2, linestyle=linestyle_usus_T0250_Cp0300
    )
ax.plot(
    data_usus_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_usus_T0250_Cp0400, linewidth=2, linestyle=linestyle_usus_T0250_Cp0400
    )
ax.plot(
    data_usus_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_usus_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_usus_T0250_Cp0500, linewidth=2, linestyle=linestyle_usus_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_usus_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0000
    )
ax.plot(
    data_usus_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0100
    )
ax.plot(
    data_usus_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0200
    )
ax.plot(
    data_usus_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0300
    )
ax.plot(
    data_usus_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0400
    )
ax.plot(
    data_usus_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_usus_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usus_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{us \rightarrow us} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usus_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usus_T0250_Cp0000, data_usus_T0250_Cp0100, data_usus_T0250_Cp0200, data_usus_T0250_Cp0300, data_usus_T0250_Cp0400, data_usus_T0250_Cp0500
del color_usus_T0250_Cp0000, color_usus_T0250_Cp0100, color_usus_T0250_Cp0200, color_usus_T0250_Cp0300, color_usus_T0250_Cp0400, color_usus_T0250_Cp0500
del linestyle_usus_T0250_Cp0000, linestyle_usus_T0250_Cp0100, linestyle_usus_T0250_Cp0200, linestyle_usus_T0250_Cp0300, linestyle_usus_T0250_Cp0400, linestyle_usus_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ubarsbar->ubarsbar
print("Building plot: cross section for ubarsbar->ubarsbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ubarsbarubarsbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ubarsbarubarsbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ubarsbarubarsbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ubarsbarubarsbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ubarsbarubarsbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUBarSBarUBarSBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ubarsbarubarsbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ubarsbarubarsbar_T0250_Cp0000 = 'black'
linestyle_ubarsbarubarsbar_T0250_Cp0000 = '-'
color_ubarsbarubarsbar_T0250_Cp0100 = 'red'
linestyle_ubarsbarubarsbar_T0250_Cp0100 = '-'
color_ubarsbarubarsbar_T0250_Cp0200 = 'blue'
linestyle_ubarsbarubarsbar_T0250_Cp0200 = '-'
color_ubarsbarubarsbar_T0250_Cp0300 = 'lime'
linestyle_ubarsbarubarsbar_T0250_Cp0300 = '-'
color_ubarsbarubarsbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ubarsbarubarsbar_T0250_Cp0400 = '-'
color_ubarsbarubarsbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ubarsbarubarsbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ubarsbarubarsbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0000
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ubarsbarubarsbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0100
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ubarsbarubarsbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0200
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ubarsbarubarsbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0300
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ubarsbarubarsbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0400
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ubarsbarubarsbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ubarsbarubarsbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_ubarsbarubarsbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0000
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0100
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0200
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0300
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0400
    )
ax.plot(
    data_ubarsbarubarsbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ubarsbarubarsbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ubarsbarubarsbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\bar{u}\bar{s} \rightarrow \bar{u}\bar{s}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 3.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add title to plot
#fig.suptitle(r'T[GeV]=0.250', fontsize=16)

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ubarsbarubarsbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ubarsbarubarsbar_T0250_Cp0000, data_ubarsbarubarsbar_T0250_Cp0100, data_ubarsbarubarsbar_T0250_Cp0200, data_ubarsbarubarsbar_T0250_Cp0300, data_ubarsbarubarsbar_T0250_Cp0400, data_ubarsbarubarsbar_T0250_Cp0500
del color_ubarsbarubarsbar_T0250_Cp0000, color_ubarsbarubarsbar_T0250_Cp0100, color_ubarsbarubarsbar_T0250_Cp0200, color_ubarsbarubarsbar_T0250_Cp0300, color_ubarsbarubarsbar_T0250_Cp0400, color_ubarsbarubarsbar_T0250_Cp0500
del linestyle_ubarsbarubarsbar_T0250_Cp0000, linestyle_ubarsbarubarsbar_T0250_Cp0100, linestyle_ubarsbarubarsbar_T0250_Cp0200, linestyle_ubarsbarubarsbar_T0250_Cp0300, linestyle_ubarsbarubarsbar_T0250_Cp0400, linestyle_ubarsbarubarsbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ss->ss
print("Building plot: cross section for ss->ss (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionSSSS_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssss_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ssss_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ssss_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ssss_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ssss_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSSS_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ssss_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ssss_T0250_Cp0000 = 'black'
linestyle_ssss_T0250_Cp0000 = '-'
color_ssss_T0250_Cp0100 = 'red'
linestyle_ssss_T0250_Cp0100 = '-'
color_ssss_T0250_Cp0200 = 'blue'
linestyle_ssss_T0250_Cp0200 = '-'
color_ssss_T0250_Cp0300 = 'lime'
linestyle_ssss_T0250_Cp0300 = '-'
color_ssss_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ssss_T0250_Cp0400 = '-'
color_ssss_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ssss_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ssss_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ssss_T0250_Cp0000, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0000
    )
ax.plot(
    data_ssss_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ssss_T0250_Cp0100, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0100
    )
ax.plot(
    data_ssss_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ssss_T0250_Cp0200, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0200
    )
ax.plot(
    data_ssss_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ssss_T0250_Cp0300, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0300
    )
ax.plot(
    data_ssss_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ssss_T0250_Cp0400, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0400
    )
ax.plot(
    data_ssss_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ssss_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ssss_T0250_Cp0500, linewidth=2, linestyle=linestyle_ssss_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ssss_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0000
    )
ax.plot(
    data_ssss_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0100
    )
ax.plot(
    data_ssss_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0200
    )
ax.plot(
    data_ssss_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0300
    )
ax.plot(
    data_ssss_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0400
    )
ax.plot(
    data_ssss_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ssss_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssss_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{ss \rightarrow ss} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ssss_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ssss_T0250_Cp0000, data_ssss_T0250_Cp0100, data_ssss_T0250_Cp0200, data_ssss_T0250_Cp0300, data_ssss_T0250_Cp0400, data_ssss_T0250_Cp0500
del color_ssss_T0250_Cp0000, color_ssss_T0250_Cp0100, color_ssss_T0250_Cp0200, color_ssss_T0250_Cp0300, color_ssss_T0250_Cp0400, color_ssss_T0250_Cp0500
del linestyle_ssss_T0250_Cp0000, linestyle_ssss_T0250_Cp0100, linestyle_ssss_T0250_Cp0200, linestyle_ssss_T0250_Cp0300, linestyle_ssss_T0250_Cp0400, linestyle_ssss_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for sbarsbar->sbarsbar
print("Building plot: cross section for sbarsbar->sbarsbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_sbarsbarsbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_sbarsbarsbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_sbarsbarsbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_sbarsbarsbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_sbarsbarsbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionSBarSBarSBarSBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_sbarsbarsbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_sbarsbarsbar_T0250_Cp0000 = 'black'
linestyle_sbarsbarsbar_T0250_Cp0000 = '-'
color_sbarsbarsbar_T0250_Cp0100 = 'red'
linestyle_sbarsbarsbar_T0250_Cp0100 = '-'
color_sbarsbarsbar_T0250_Cp0200 = 'blue'
linestyle_sbarsbarsbar_T0250_Cp0200 = '-'
color_sbarsbarsbar_T0250_Cp0300 = 'lime'
linestyle_sbarsbarsbar_T0250_Cp0300 = '-'
color_sbarsbarsbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_sbarsbarsbar_T0250_Cp0400 = '-'
color_sbarsbarsbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_sbarsbarsbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_sbarsbarsbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_sbarsbarsbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0000
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_sbarsbarsbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0100
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_sbarsbarsbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0200
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_sbarsbarsbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0300
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_sbarsbarsbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0400
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_sbarsbarsbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_sbarsbarsbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_sbarsbarsbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_sbarsbarsbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0000
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0100
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0200
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0300
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0400
    )
ax.plot(
    data_sbarsbarsbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_sbarsbarsbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_sbarsbarsbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\bar{s}\bar{s} \rightarrow \bar{s}\bar{s}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 5.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=6, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_sbarsbarsbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_sbarsbarsbar_T0250_Cp0000, data_sbarsbarsbar_T0250_Cp0100, data_sbarsbarsbar_T0250_Cp0200, data_sbarsbarsbar_T0250_Cp0300, data_sbarsbarsbar_T0250_Cp0400, data_sbarsbarsbar_T0250_Cp0500
del color_sbarsbarsbar_T0250_Cp0000, color_sbarsbarsbar_T0250_Cp0100, color_sbarsbarsbar_T0250_Cp0200, color_sbarsbarsbar_T0250_Cp0300, color_sbarsbarsbar_T0250_Cp0400, color_sbarsbarsbar_T0250_Cp0500
del linestyle_sbarsbarsbar_T0250_Cp0000, linestyle_sbarsbarsbar_T0250_Cp0100, linestyle_sbarsbarsbar_T0250_Cp0200, linestyle_sbarsbarsbar_T0250_Cp0300, linestyle_sbarsbarsbar_T0250_Cp0400, linestyle_sbarsbarsbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for uubar->uubar
print("Building plot: cross section for uubar->uubar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUUBarUUBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubaruubar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarUUBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_uubaruubar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarUUBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_uubaruubar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarUUBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_uubaruubar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarUUBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_uubaruubar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarUUBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_uubaruubar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uubaruubar_T0250_Cp0000 = 'black'
linestyle_uubaruubar_T0250_Cp0000 = '-'
color_uubaruubar_T0250_Cp0100 = 'red'
linestyle_uubaruubar_T0250_Cp0100 = '-'
color_uubaruubar_T0250_Cp0200 = 'blue'
linestyle_uubaruubar_T0250_Cp0200 = '-'
color_uubaruubar_T0250_Cp0300 = 'lime'
linestyle_uubaruubar_T0250_Cp0300 = '-'
color_uubaruubar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_uubaruubar_T0250_Cp0400 = '-'
color_uubaruubar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_uubaruubar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_uubaruubar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_uubaruubar_T0250_Cp0000, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0000
    )
ax.plot(
    data_uubaruubar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_uubaruubar_T0250_Cp0100, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0100
    )
ax.plot(
    data_uubaruubar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_uubaruubar_T0250_Cp0200, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0200
    )
ax.plot(
    data_uubaruubar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_uubaruubar_T0250_Cp0300, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0300
    )
ax.plot(
    data_uubaruubar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_uubaruubar_T0250_Cp0400, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0400
    )
ax.plot(
    data_uubaruubar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_uubaruubar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_uubaruubar_T0250_Cp0500, linewidth=2, linestyle=linestyle_uubaruubar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_uubaruubar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0000
    )
ax.plot(
    data_uubaruubar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0100
    )
ax.plot(
    data_uubaruubar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0200
    )
ax.plot(
    data_uubaruubar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0300
    )
ax.plot(
    data_uubaruubar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0400
    )
ax.plot(
    data_uubaruubar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_uubaruubar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubaruubar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{u\bar{u} \rightarrow u\bar{u}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 40.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uubaruubar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uubaruubar_T0250_Cp0000, data_uubaruubar_T0250_Cp0100, data_uubaruubar_T0250_Cp0200, data_uubaruubar_T0250_Cp0300, data_uubaruubar_T0250_Cp0400, data_uubaruubar_T0250_Cp0500
del color_uubaruubar_T0250_Cp0000, color_uubaruubar_T0250_Cp0100, color_uubaruubar_T0250_Cp0200, color_uubaruubar_T0250_Cp0300, color_uubaruubar_T0250_Cp0400, color_uubaruubar_T0250_Cp0500
del linestyle_uubaruubar_T0250_Cp0000, linestyle_uubaruubar_T0250_Cp0100, linestyle_uubaruubar_T0250_Cp0200, linestyle_uubaruubar_T0250_Cp0300, linestyle_uubaruubar_T0250_Cp0400, linestyle_uubaruubar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for udbar->udbar
print("Building plot: cross section for udbar->udbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUDBarUDBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_udbarudbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDBarUDBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_udbarudbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDBarUDBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_udbarudbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDBarUDBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_udbarudbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDBarUDBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_udbarudbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUDBarUDBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_udbarudbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_udbarudbar_T0250_Cp0000 = 'black'
linestyle_udbarudbar_T0250_Cp0000 = '-'
color_udbarudbar_T0250_Cp0100 = 'red'
linestyle_udbarudbar_T0250_Cp0100 = '-'
color_udbarudbar_T0250_Cp0200 = 'blue'
linestyle_udbarudbar_T0250_Cp0200 = '-'
color_udbarudbar_T0250_Cp0300 = 'lime'
linestyle_udbarudbar_T0250_Cp0300 = '-'
color_udbarudbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_udbarudbar_T0250_Cp0400 = '-'
color_udbarudbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_udbarudbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_udbarudbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_udbarudbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0000
    )
ax.plot(
    data_udbarudbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_udbarudbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0100
    )
ax.plot(
    data_udbarudbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_udbarudbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0200
    )
ax.plot(
    data_udbarudbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_udbarudbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0300
    )
ax.plot(
    data_udbarudbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_udbarudbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0400
    )
ax.plot(
    data_udbarudbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_udbarudbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_udbarudbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_udbarudbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_udbarudbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0000
    )
ax.plot(
    data_udbarudbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0100
    )
ax.plot(
    data_udbarudbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0200
    )
ax.plot(
    data_udbarudbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0300
    )
ax.plot(
    data_udbarudbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0400
    )
ax.plot(
    data_udbarudbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_udbarudbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_udbarudbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{u\bar{d} \rightarrow u\bar{d}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 40.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_udbarudbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_udbarudbar_T0250_Cp0000, data_udbarudbar_T0250_Cp0100, data_udbarudbar_T0250_Cp0200, data_udbarudbar_T0250_Cp0300, data_udbarudbar_T0250_Cp0400, data_udbarudbar_T0250_Cp0500
del color_udbarudbar_T0250_Cp0000, color_udbarudbar_T0250_Cp0100, color_udbarudbar_T0250_Cp0200, color_udbarudbar_T0250_Cp0300, color_udbarudbar_T0250_Cp0400, color_udbarudbar_T0250_Cp0500
del linestyle_udbarudbar_T0250_Cp0000, linestyle_udbarudbar_T0250_Cp0100, linestyle_udbarudbar_T0250_Cp0200, linestyle_udbarudbar_T0250_Cp0300, linestyle_udbarudbar_T0250_Cp0400, linestyle_udbarudbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for uubar->ddbar
print("Building plot: cross section for uubar->ddbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUUBarDDBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarddbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarDDBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_uubarddbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarDDBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_uubarddbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarDDBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_uubarddbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarDDBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_uubarddbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarDDBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_uubarddbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uubarddbar_T0250_Cp0000 = 'black'
linestyle_uubarddbar_T0250_Cp0000 = '-'
color_uubarddbar_T0250_Cp0100 = 'red'
linestyle_uubarddbar_T0250_Cp0100 = '-'
color_uubarddbar_T0250_Cp0200 = 'blue'
linestyle_uubarddbar_T0250_Cp0200 = '-'
color_uubarddbar_T0250_Cp0300 = 'lime'
linestyle_uubarddbar_T0250_Cp0300 = '-'
color_uubarddbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_uubarddbar_T0250_Cp0400 = '-'
color_uubarddbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_uubarddbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_uubarddbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_uubarddbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0000
    )
ax.plot(
    data_uubarddbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_uubarddbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0100
    )
ax.plot(
    data_uubarddbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_uubarddbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0200
    )
ax.plot(
    data_uubarddbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_uubarddbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0300
    )
ax.plot(
    data_uubarddbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_uubarddbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0400
    )
ax.plot(
    data_uubarddbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_uubarddbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_uubarddbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_uubarddbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_uubarddbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0000
    )
ax.plot(
    data_uubarddbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0100
    )
ax.plot(
    data_uubarddbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0200
    )
ax.plot(
    data_uubarddbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0300
    )
ax.plot(
    data_uubarddbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0400
    )
ax.plot(
    data_uubarddbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_uubarddbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarddbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{u\bar{u} \rightarrow d\bar{d}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper right', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 40.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uubarddbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uubarddbar_T0250_Cp0000, data_uubarddbar_T0250_Cp0100, data_uubarddbar_T0250_Cp0200, data_uubarddbar_T0250_Cp0300, data_uubarddbar_T0250_Cp0400, data_uubarddbar_T0250_Cp0500
del color_uubarddbar_T0250_Cp0000, color_uubarddbar_T0250_Cp0100, color_uubarddbar_T0250_Cp0200, color_uubarddbar_T0250_Cp0300, color_uubarddbar_T0250_Cp0400, color_uubarddbar_T0250_Cp0500
del linestyle_uubarddbar_T0250_Cp0000, linestyle_uubarddbar_T0250_Cp0100, linestyle_uubarddbar_T0250_Cp0200, linestyle_uubarddbar_T0250_Cp0300, linestyle_uubarddbar_T0250_Cp0400, linestyle_uubarddbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for uubar->ssbar
print("Building plot: cross section for uubar->ssbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUUBarSSBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_uubarssbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarSSBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_uubarssbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarSSBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_uubarssbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarSSBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_uubarssbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarSSBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_uubarssbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUUBarSSBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_uubarssbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_uubarssbar_T0250_Cp0000 = 'black'
linestyle_uubarssbar_T0250_Cp0000 = '-'
color_uubarssbar_T0250_Cp0100 = 'red'
linestyle_uubarssbar_T0250_Cp0100 = '-'
color_uubarssbar_T0250_Cp0200 = 'blue'
linestyle_uubarssbar_T0250_Cp0200 = '-'
color_uubarssbar_T0250_Cp0300 = 'lime'
linestyle_uubarssbar_T0250_Cp0300 = '-'
color_uubarssbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_uubarssbar_T0250_Cp0400 = '-'
color_uubarssbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_uubarssbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_uubarssbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_uubarssbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0000
    )
ax.plot(
    data_uubarssbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_uubarssbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0100
    )
ax.plot(
    data_uubarssbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_uubarssbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0200
    )
ax.plot(
    data_uubarssbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_uubarssbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0300
    )
ax.plot(
    data_uubarssbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_uubarssbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0400
    )
ax.plot(
    data_uubarssbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_uubarssbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_uubarssbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_uubarssbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_uubarssbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0000
    )
ax.plot(
    data_uubarssbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0100
    )
ax.plot(
    data_uubarssbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0200
    )
ax.plot(
    data_uubarssbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0300
    )
ax.plot(
    data_uubarssbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0400
    )
ax.plot(
    data_uubarssbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_uubarssbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_uubarssbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{u\bar{u} \rightarrow s\bar{s}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 2.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_uubarssbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_uubarssbar_T0250_Cp0000, data_uubarssbar_T0250_Cp0100, data_uubarssbar_T0250_Cp0200, data_uubarssbar_T0250_Cp0300, data_uubarssbar_T0250_Cp0400, data_uubarssbar_T0250_Cp0500
del color_uubarssbar_T0250_Cp0000, color_uubarssbar_T0250_Cp0100, color_uubarssbar_T0250_Cp0200, color_uubarssbar_T0250_Cp0300, color_uubarssbar_T0250_Cp0400, color_uubarssbar_T0250_Cp0500
del linestyle_uubarssbar_T0250_Cp0000, linestyle_uubarssbar_T0250_Cp0100, linestyle_uubarssbar_T0250_Cp0200, linestyle_uubarssbar_T0250_Cp0300, linestyle_uubarssbar_T0250_Cp0400, linestyle_uubarssbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ssbar->uubar
print("Building plot: cross section for ssbar->uubar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionSSBarUUBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbaruubar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarUUBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ssbaruubar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarUUBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ssbaruubar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarUUBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ssbaruubar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarUUBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ssbaruubar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarUUBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ssbaruubar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ssbaruubar_T0250_Cp0000 = 'black'
linestyle_ssbaruubar_T0250_Cp0000 = '-'
color_ssbaruubar_T0250_Cp0100 = 'red'
linestyle_ssbaruubar_T0250_Cp0100 = '-'
color_ssbaruubar_T0250_Cp0200 = 'blue'
linestyle_ssbaruubar_T0250_Cp0200 = '-'
color_ssbaruubar_T0250_Cp0300 = 'lime'
linestyle_ssbaruubar_T0250_Cp0300 = '-'
color_ssbaruubar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ssbaruubar_T0250_Cp0400 = '-'
color_ssbaruubar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ssbaruubar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ssbaruubar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ssbaruubar_T0250_Cp0000, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0000
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ssbaruubar_T0250_Cp0100, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0100
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ssbaruubar_T0250_Cp0200, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0200
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ssbaruubar_T0250_Cp0300, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0300
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ssbaruubar_T0250_Cp0400, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0400
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ssbaruubar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ssbaruubar_T0250_Cp0500, linewidth=2, linestyle=linestyle_ssbaruubar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ssbaruubar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0000
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0100
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0200
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0300
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0400
    )
ax.plot(
    data_ssbaruubar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ssbaruubar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbaruubar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{s\bar{s} \rightarrow u\bar{u}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 20.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ssbaruubar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ssbaruubar_T0250_Cp0000, data_ssbaruubar_T0250_Cp0100, data_ssbaruubar_T0250_Cp0200, data_ssbaruubar_T0250_Cp0300, data_ssbaruubar_T0250_Cp0400, data_ssbaruubar_T0250_Cp0500
del color_ssbaruubar_T0250_Cp0000, color_ssbaruubar_T0250_Cp0100, color_ssbaruubar_T0250_Cp0200, color_ssbaruubar_T0250_Cp0300, color_ssbaruubar_T0250_Cp0400, color_ssbaruubar_T0250_Cp0500
del linestyle_ssbaruubar_T0250_Cp0000, linestyle_ssbaruubar_T0250_Cp0100, linestyle_ssbaruubar_T0250_Cp0200, linestyle_ssbaruubar_T0250_Cp0300, linestyle_ssbaruubar_T0250_Cp0400, linestyle_ssbaruubar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for ssbar->ssbar
print("Building plot: cross section for ssbar->ssbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionSSBarSSBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_ssbarssbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarSSBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_ssbarssbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarSSBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_ssbarssbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarSSBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_ssbarssbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarSSBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_ssbarssbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionSSBarSSBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_ssbarssbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_ssbarssbar_T0250_Cp0000 = 'black'
linestyle_ssbarssbar_T0250_Cp0000 = '-'
color_ssbarssbar_T0250_Cp0100 = 'red'
linestyle_ssbarssbar_T0250_Cp0100 = '-'
color_ssbarssbar_T0250_Cp0200 = 'blue'
linestyle_ssbarssbar_T0250_Cp0200 = '-'
color_ssbarssbar_T0250_Cp0300 = 'lime'
linestyle_ssbarssbar_T0250_Cp0300 = '-'
color_ssbarssbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_ssbarssbar_T0250_Cp0400 = '-'
color_ssbarssbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_ssbarssbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_ssbarssbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_ssbarssbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0000
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_ssbarssbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0100
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_ssbarssbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0200
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_ssbarssbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0300
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_ssbarssbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0400
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_ssbarssbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_ssbarssbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_ssbarssbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_ssbarssbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0000
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0100
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0200
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0300
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0400
    )
ax.plot(
    data_ssbarssbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_ssbarssbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_ssbarssbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{s\bar{s} \rightarrow s\bar{s}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 30.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_ssbarssbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_ssbarssbar_T0250_Cp0000, data_ssbarssbar_T0250_Cp0100, data_ssbarssbar_T0250_Cp0200, data_ssbarssbar_T0250_Cp0300, data_ssbarssbar_T0250_Cp0400, data_ssbarssbar_T0250_Cp0500
del color_ssbarssbar_T0250_Cp0000, color_ssbarssbar_T0250_Cp0100, color_ssbarssbar_T0250_Cp0200, color_ssbarssbar_T0250_Cp0300, color_ssbarssbar_T0250_Cp0400, color_ssbarssbar_T0250_Cp0500
del linestyle_ssbarssbar_T0250_Cp0000, linestyle_ssbarssbar_T0250_Cp0100, linestyle_ssbarssbar_T0250_Cp0200, linestyle_ssbarssbar_T0250_Cp0300, linestyle_ssbarssbar_T0250_Cp0400, linestyle_ssbarssbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for usbar->usbar
print("Building plot: cross section for usbar->usbar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionUSBarUSBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_usbarusbar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSBarUSBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_usbarusbar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSBarUSBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_usbarusbar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSBarUSBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_usbarusbar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSBarUSBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_usbarusbar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionUSBarUSBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_usbarusbar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_usbarusbar_T0250_Cp0000 = 'black'
linestyle_usbarusbar_T0250_Cp0000 = '-'
color_usbarusbar_T0250_Cp0100 = 'red'
linestyle_usbarusbar_T0250_Cp0100 = '-'
color_usbarusbar_T0250_Cp0200 = 'blue'
linestyle_usbarusbar_T0250_Cp0200 = '-'
color_usbarusbar_T0250_Cp0300 = 'lime'
linestyle_usbarusbar_T0250_Cp0300 = '-'
color_usbarusbar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_usbarusbar_T0250_Cp0400 = '-'
color_usbarusbar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_usbarusbar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_usbarusbar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_usbarusbar_T0250_Cp0000, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0000
    )
ax.plot(
    data_usbarusbar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_usbarusbar_T0250_Cp0100, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0100
    )
ax.plot(
    data_usbarusbar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_usbarusbar_T0250_Cp0200, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0200
    )
ax.plot(
    data_usbarusbar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_usbarusbar_T0250_Cp0300, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0300
    )
ax.plot(
    data_usbarusbar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_usbarusbar_T0250_Cp0400, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0400
    )
ax.plot(
    data_usbarusbar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_usbarusbar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_usbarusbar_T0250_Cp0500, linewidth=2, linestyle=linestyle_usbarusbar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_usbarusbar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0000
    )
ax.plot(
    data_usbarusbar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0100
    )
ax.plot(
    data_usbarusbar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0200
    )
ax.plot(
    data_usbarusbar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0300
    )
ax.plot(
    data_usbarusbar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0400
    )
ax.plot(
    data_usbarusbar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_usbarusbar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_usbarusbar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{u\bar{s} \rightarrow u\bar{s}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 30.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_usbarusbar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_usbarusbar_T0250_Cp0000, data_usbarusbar_T0250_Cp0100, data_usbarusbar_T0250_Cp0200, data_usbarusbar_T0250_Cp0300, data_usbarusbar_T0250_Cp0400, data_usbarusbar_T0250_Cp0500
del color_usbarusbar_T0250_Cp0000, color_usbarusbar_T0250_Cp0100, color_usbarusbar_T0250_Cp0200, color_usbarusbar_T0250_Cp0300, color_usbarusbar_T0250_Cp0400, color_usbarusbar_T0250_Cp0500
del linestyle_usbarusbar_T0250_Cp0000, linestyle_usbarusbar_T0250_Cp0100, linestyle_usbarusbar_T0250_Cp0200, linestyle_usbarusbar_T0250_Cp0300, linestyle_usbarusbar_T0250_Cp0400, linestyle_usbarusbar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# Cross section for subar->subar
print("Building plot: cross section for subar->subar (T[GeV]=0.250, finite chemical potential)")

# Load data from the file
filename = "crossSectionSUBarSUBar_T0.250000_CPU0.000000_CPD0.000000_CPS0.000000.dat"
data_subarsubar_T0250_Cp0000 = CrossSectionData(data_folder + filename)

filename = "crossSectionSUBarSUBar_T0.250000_CPU0.100000_CPD0.100000_CPS0.100000.dat"
data_subarsubar_T0250_Cp0100 = CrossSectionData(data_folder + filename)

filename = "crossSectionSUBarSUBar_T0.250000_CPU0.200000_CPD0.200000_CPS0.200000.dat"
data_subarsubar_T0250_Cp0200 = CrossSectionData(data_folder + filename)

filename = "crossSectionSUBarSUBar_T0.250000_CPU0.300000_CPD0.300000_CPS0.300000.dat"
data_subarsubar_T0250_Cp0300 = CrossSectionData(data_folder + filename)

filename = "crossSectionSUBarSUBar_T0.250000_CPU0.400000_CPD0.400000_CPS0.400000.dat"
data_subarsubar_T0250_Cp0400 = CrossSectionData(data_folder + filename)

filename = "crossSectionSUBarSUBar_T0.250000_CPU0.500000_CPD0.500000_CPS0.500000.dat"
data_subarsubar_T0250_Cp0500 = CrossSectionData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_subarsubar_T0250_Cp0000 = 'black'
linestyle_subarsubar_T0250_Cp0000 = '-'
color_subarsubar_T0250_Cp0100 = 'red'
linestyle_subarsubar_T0250_Cp0100 = '-'
color_subarsubar_T0250_Cp0200 = 'blue'
linestyle_subarsubar_T0250_Cp0200 = '-'
color_subarsubar_T0250_Cp0300 = 'lime'
linestyle_subarsubar_T0250_Cp0300 = '-'
color_subarsubar_T0250_Cp0400 = (0, 0.9, 0.9)  # Light cyan/turquoise
linestyle_subarsubar_T0250_Cp0400 = '-'
color_subarsubar_T0250_Cp0500 = (1, 0.1, 0.9)   # Magenta
linestyle_subarsubar_T0250_Cp0500 = '-'

# Plot cross sections
ax.plot(
    data_subarsubar_T0250_Cp0000.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0000.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.0$', color=color_subarsubar_T0250_Cp0000, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0000
    )
ax.plot(
    data_subarsubar_T0250_Cp0100.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0100.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.1$', color=color_subarsubar_T0250_Cp0100, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0100
    )
ax.plot(
    data_subarsubar_T0250_Cp0200.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0200.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.2$', color=color_subarsubar_T0250_Cp0200, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0200
    )
ax.plot(
    data_subarsubar_T0250_Cp0300.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0300.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.3$', color=color_subarsubar_T0250_Cp0300, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0300
    )
ax.plot(
    data_subarsubar_T0250_Cp0400.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0400.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.4$', color=color_subarsubar_T0250_Cp0400, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0400
    )
ax.plot(
    data_subarsubar_T0250_Cp0500.get_sqrt_center_of_mass_energy(), 
    data_subarsubar_T0250_Cp0500.get_cross_section(), 
    label=r'$\mu_q [\mathrm{GeV}]=0.5$', color=color_subarsubar_T0250_Cp0500, linewidth=2, linestyle=linestyle_subarsubar_T0250_Cp0500
    )

# Plot minimum center of mass energy point
ax.plot(
    data_subarsubar_T0250_Cp0000.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0000.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0000
    )
ax.plot(
    data_subarsubar_T0250_Cp0100.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0100.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0100
    )
ax.plot(
    data_subarsubar_T0250_Cp0200.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0200.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0200
    )
ax.plot(
    data_subarsubar_T0250_Cp0300.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0300.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0300
    )
ax.plot(
    data_subarsubar_T0250_Cp0400.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0400.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0400
    )
ax.plot(
    data_subarsubar_T0250_Cp0500.get_sqrt_center_of_mass_energy()[0], 
    data_subarsubar_T0250_Cp0500.get_cross_section()[0], 
    marker='o', markersize=6, color=color_subarsubar_T0250_Cp0500
    )

# Axes labels
ax.set_xlabel(r'$\sqrt{s} \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{s\bar{u} \rightarrow s\bar{u}} \, [\mathrm{mb}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False, title=r'T[GeV]=0.250', title_fontsize=14)

# Configure axes using the helper function
xmin= 0; xmax = 1.2; ymin = 0.0; ymax = 30.0
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=7, y_num_ticks=7, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "cross_section_subarsubar_T0250_CP.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_subarsubar_T0250_Cp0000, data_subarsubar_T0250_Cp0100, data_subarsubar_T0250_Cp0200, data_subarsubar_T0250_Cp0300, data_subarsubar_T0250_Cp0400, data_subarsubar_T0250_Cp0500
del color_subarsubar_T0250_Cp0000, color_subarsubar_T0250_Cp0100, color_subarsubar_T0250_Cp0200, color_subarsubar_T0250_Cp0300, color_subarsubar_T0250_Cp0400, color_subarsubar_T0250_Cp0500
del linestyle_subarsubar_T0250_Cp0000, linestyle_subarsubar_T0250_Cp0100, linestyle_subarsubar_T0250_Cp0200, linestyle_subarsubar_T0250_Cp0300, linestyle_subarsubar_T0250_Cp0400, linestyle_subarsubar_T0250_Cp0500
del fig, ax, xmin, xmax, ymin, ymax

