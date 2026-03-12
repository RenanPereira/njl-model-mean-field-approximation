import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import *
from common_utils.first_order_line_data import *


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
data_folder = "su3_3d_cutoff_phase_diagram/data/"
plots_folder = "su3_3d_cutoff_phase_diagram/plots/"


####################################################################################################
# First order line set A
print("Building plot: first order line (chemical potential versus temperature) for set A")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setA.dat"
data_setA = FirstOrderLineData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setA = 'black'

# Plot first order line
ax.plot(data_setA.get_quark_chem_pot(), data_setA.get_temperature(), label=r'Set A', color=color_setA, linewidth=2, linestyle='-')

# Plot Critical End point
ax.plot(data_setA.get_quark_chem_pot()[-1], data_setA.get_temperature()[-1], marker='o', markersize=10, color=color_setA)

# Axes labels
ax.set_xlabel(r'$\mu_q \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.4; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "first_order_setA_muq_vs_temp.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setA
del color_setA
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# First order line set B
print("Building plot: first order line (chemical potential versus temperature) for set B")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setB.dat"
data_setB = FirstOrderLineData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setB = 'red'

# Plot first order line
ax.plot(data_setB.get_quark_chem_pot(), data_setB.get_temperature(), label=r'Set B', color=color_setB, linewidth=2, linestyle='-')

# Plot Critical End point
ax.plot(data_setB.get_quark_chem_pot()[-1], data_setB.get_temperature()[-1], marker='o', markersize=10, color=color_setB)

# Axes labels
ax.set_xlabel(r'$\mu_q \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.4; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "first_order_setB_muq_vs_temp.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setB
del color_setB
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# First order line set C
print("Building plot: first order line (chemical potential versus temperature) for set C")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setC.dat"
data_setC = FirstOrderLineData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setC = 'blue'

# Plot first order line
ax.plot(data_setC.get_quark_chem_pot(), data_setC.get_temperature(), label=r'Set C', color=color_setC, linewidth=2, linestyle='-')

# Plot Critical End point
ax.plot(data_setC.get_quark_chem_pot()[-1], data_setC.get_temperature()[-1], marker='o', markersize=10, color=color_setC)

# Axes labels
ax.set_xlabel(r'$\mu_q \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.4; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "first_order_setC_muq_vs_temp.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setC
del color_setC
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# First order line all sets
print("Building plot: first order line (chemical potential versus temperature) for sets A, B and C")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setA.dat"
data_setA = FirstOrderLineData(data_folder + filename)

filename = "SU3NJL3DCutoffFirstOrderLine_setB.dat"
data_setB = FirstOrderLineData(data_folder + filename)

filename = "SU3NJL3DCutoffFirstOrderLine_setC.dat"
data_setC = FirstOrderLineData(data_folder + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setA = 'black'
color_setB = 'red'
color_setC = 'blue'

# Plot first order line
ax.plot(data_setA.get_quark_chem_pot(), data_setA.get_temperature(), label=r'Set A', color=color_setA, linewidth=2, linestyle='-')
ax.plot(data_setB.get_quark_chem_pot(), data_setB.get_temperature(), label=r'Set B', color=color_setB, linewidth=2, linestyle='-')
ax.plot(data_setC.get_quark_chem_pot(), data_setC.get_temperature(), label=r'Set C', color=color_setC, linewidth=2, linestyle='-')

# Plot Critical End point
ax.plot(data_setA.get_quark_chem_pot()[-1], data_setA.get_temperature()[-1], marker='o', markersize=10, color=color_setA)
ax.plot(data_setB.get_quark_chem_pot()[-1], data_setB.get_temperature()[-1], marker='o', markersize=10, color=color_setB)
ax.plot(data_setC.get_quark_chem_pot()[-1], data_setC.get_temperature()[-1], marker='o', markersize=10, color=color_setC)

# Axes labels
ax.set_xlabel(r'$\mu_q \, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.4; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=5, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "first_order_allSets_muq_vs_temp.png"
plt.savefig(plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setA, data_setB, data_setC
del color_setA, color_setB, color_setC
del fig, ax, xmin, xmax, ymin, ymax
