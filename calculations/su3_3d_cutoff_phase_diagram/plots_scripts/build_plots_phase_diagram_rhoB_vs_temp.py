import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import sys

# Get the parent directory of the current script (two levels up) and add the parent directory to sys.path
# This is necessary in order to import the PhaseTransitionAnalyzerAtFixedTemperature class
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'common_utils'))
sys.path.append(parent_dir)

from plot_helper import *
from first_order_line_data import *


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6


####################################################################################################
# Phase Diagram set A
print("Building plot: first order line (baryon density versus temperature) for set A")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setA.dat"
data_setA = FirstOrderLineData( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setA = 'black'

# Plot first order line
ax.plot(data_setA.get_baryon_density_broken(), data_setA.get_temperature(), label=r'Set A', color=color_setA, linewidth=2, linestyle='-')
ax.plot(data_setA.get_baryon_density_restored(), data_setA.get_temperature(), color=color_setA, linewidth=2, linestyle='-')

# Plot Critical End point
cep_rhoB = 0.5*(data_setA.get_baryon_density_broken()[-1] + data_setA.get_baryon_density_restored()[-1])
cep_temp = data_setA.get_temperature()[-1]
ax.plot(cep_rhoB, cep_temp, marker='o', markersize=10, color=color_setA)

# Axes labels
ax.set_xlabel(r'$\rho_{\mathrm{B}} \, [\mathrm{fm}^{-3}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.45; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=4, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "first_order_setA_rhoB_vs_temp.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setA
del color_setA
del cep_rhoB, cep_temp 
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# First order line set B
print("Building plot: first order line (baryon density versus temperature) for set B")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setB.dat"
data_setB = FirstOrderLineData( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setB = 'red'

# Plot first order line
ax.plot(data_setB.get_baryon_density_broken(), data_setB.get_temperature(), label=r'Set B', color=color_setB, linewidth=2, linestyle='-')
ax.plot(data_setB.get_baryon_density_restored(), data_setB.get_temperature(), color=color_setB, linewidth=2, linestyle='-')

# Plot Critical End point
cep_rhoB = 0.5*(data_setB.get_baryon_density_broken()[-1] + data_setB.get_baryon_density_restored()[-1])
cep_temp = data_setB.get_temperature()[-1]
ax.plot(cep_rhoB, cep_temp, marker='o', markersize=10, color=color_setB)

# Axes labels
ax.set_xlabel(r'$\rho_{\mathrm{B}} \, [\mathrm{fm}^{-3}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.45; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=4, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "first_order_setB_rhoB_vs_temp.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setB
del color_setB
del cep_rhoB, cep_temp 
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# First order line set C
print("Building plot: first order line (baryon density versus temperature) for set C")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setC.dat"
data_setC = FirstOrderLineData( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setC = 'blue'

# Plot first order line
ax.plot(data_setC.get_baryon_density_broken(), data_setC.get_temperature(), label=r'Set C', color=color_setC, linewidth=2, linestyle='-')
ax.plot(data_setC.get_baryon_density_restored(), data_setC.get_temperature(), color=color_setC, linewidth=2, linestyle='-')

# Plot Critical End point
cep_rhoB = 0.5*(data_setC.get_baryon_density_broken()[-1] + data_setC.get_baryon_density_restored()[-1])
cep_temp = data_setC.get_temperature()[-1]
ax.plot(cep_rhoB, cep_temp, marker='o', markersize=10, color=color_setC)

# Axes labels
ax.set_xlabel(r'$\rho_{\mathrm{B}} \, [\mathrm{fm}^{-3}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.45; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=4, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

# Replace .dat with .png
plotname = "first_order_setC_rhoB_vs_temp.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setC
del color_setC
del cep_rhoB, cep_temp 
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# First order line - all sets
print("Building plot: first order line (chemical potential versus temperature) for sets A, B and C")

# Load data from the file
filename = "SU3NJL3DCutoffFirstOrderLine_setA.dat"
data_setA = FirstOrderLineData( "../data/" + filename)

filename = "SU3NJL3DCutoffFirstOrderLine_setB.dat"
data_setB = FirstOrderLineData( "../data/" + filename)

filename = "SU3NJL3DCutoffFirstOrderLine_setC.dat"
data_setC = FirstOrderLineData( "../data/" + filename)

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_setA = 'black'
color_setB = 'red'
color_setC = 'blue'

# Plot first order line
ax.plot(data_setA.get_baryon_density_broken(), data_setA.get_temperature(), label=r'Set A', color=color_setA, linewidth=2, linestyle='-')
ax.plot(data_setA.get_baryon_density_restored(), data_setA.get_temperature(), color=color_setA, linewidth=2, linestyle='-')

ax.plot(data_setB.get_baryon_density_broken(), data_setB.get_temperature(), label=r'Set B', color=color_setB, linewidth=2, linestyle='-')
ax.plot(data_setB.get_baryon_density_restored(), data_setB.get_temperature(), color=color_setB, linewidth=2, linestyle='-')

ax.plot(data_setC.get_baryon_density_broken(), data_setC.get_temperature(), label=r'Set C', color=color_setC, linewidth=2, linestyle='-')
ax.plot(data_setC.get_baryon_density_restored(), data_setC.get_temperature(), color=color_setC, linewidth=2, linestyle='-')

# Plot Critical End point
cep_rhoB_setA = 0.5*(data_setA.get_baryon_density_broken()[-1] + data_setA.get_baryon_density_restored()[-1])
cep_temp_setA = data_setA.get_temperature()[-1]
ax.plot(cep_rhoB_setA, cep_temp_setA, marker='o', markersize=10, color=color_setA)

cep_rhoB_setB = 0.5*(data_setB.get_baryon_density_broken()[-1] + data_setB.get_baryon_density_restored()[-1])
cep_temp_setB = data_setB.get_temperature()[-1]
ax.plot(cep_rhoB_setB, cep_temp_setB, marker='o', markersize=10, color=color_setB)

cep_rhoB_setC = 0.5*(data_setC.get_baryon_density_broken()[-1] + data_setC.get_baryon_density_restored()[-1])
cep_temp_setC = data_setC.get_temperature()[-1]
ax.plot(cep_rhoB_setC, cep_temp_setC, marker='o', markersize=10, color=color_setC)

# Axes labels
ax.set_xlabel(r'$\rho_{\mathrm{B}} \, [\mathrm{fm}^{-3}]$', fontsize=20)
ax.set_ylabel(r'$\mathrm{T} \, [\mathrm{GeV}]$', fontsize=20)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='upper left', fontsize=14, frameon=False) # Show legend

# Configure axes using the helper function
xmin= 0; xmax = 0.45; ymin = 0.0; ymax = 0.2
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks=4, y_num_ticks=5, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Automatically adjust layout
fig.tight_layout()

plotname = "first_order_allSets_rhoB_vs_temp.png"
plt.savefig("../plots/" + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del filename, plotname
del data_setA, data_setB, data_setC
del color_setA, color_setB, color_setC
del cep_rhoB_setA, cep_temp_setA, cep_rhoB_setB, cep_temp_setB, cep_rhoB_setC, cep_temp_setC
del fig, ax, xmin, xmax, ymin, ymax
