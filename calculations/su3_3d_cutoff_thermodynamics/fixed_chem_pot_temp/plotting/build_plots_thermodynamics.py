import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"


####################################################################################################
# set A
print("Building plot: entropy density versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color = 'black'
linestyle = '-'

ax.plot(
    data.get_temperature(), 
    data.get_entropy_density(), 
    label=r'Set A', 
    color=color, 
    linewidth=2, 
    linestyle=linestyle
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.300
ymin = 0.0
ymax = 0.4
x_num_ticks = 4
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'entropy_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color, linestyle
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# set A
print("Building plot: pressure versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color = 'black'
linestyle = '-'

ax.plot(
    data.get_temperature(), 
    data.get_pressure(), 
    label=r'Set A', 
    color=color, 
    linewidth=2, 
    linestyle=linestyle
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$P\, [\mathrm{GeV}^4]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.500
ymin = 0.0
ymax = 0.2
x_num_ticks = 7
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'pressure_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color, linestyle
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# set A
print("Building plot: entropy density, dPdT versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_entropy = 'black'
linestyle_entropy = '-'

color_dPdT = 'red'
linestyle_dPdT = '--'

dPdT = np.gradient(data.get_pressure(), data.get_temperature())

ax.plot(
    data.get_temperature(), 
    data.get_entropy_density(), 
    label=r'$s$', 
    color=color_entropy, 
    linewidth=2, 
    linestyle=linestyle_entropy
)

ax.plot(
    data.get_temperature(), 
    dPdT, 
    label=r'$ ({\partial P}/{\partial T})|_{\mu}$', 
    color=color_dPdT, 
    linewidth=2, 
    linestyle=linestyle_dPdT
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.300
ymin = 0.0
ymax = 0.4
x_num_ticks = 4
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'entropy_dPdT_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color_entropy, linestyle_entropy
del color_dPdT, linestyle_dPdT
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# set A
print("Building plot: pressure vs energy density.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color = 'black'
linestyle = '-'

ax.plot(
    data.get_energy_density(), 
    data.get_pressure(), 
    label=r'Set A', 
    color=color, 
    linewidth=2, 
    linestyle=linestyle
)

# # Axes labels
ax.set_xlabel(r'$\epsilon \, [\mathrm{GeV}^4]$', fontsize=20)
ax.set_ylabel(r'$P \, [\mathrm{GeV}^4]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.6
ymin = 0.0
ymax = 0.2
x_num_ticks = 5
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'pressure_vs_energy_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color, linestyle
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# set A
print("Building plot: energy density versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color = 'black'
linestyle = '-'

ax.plot(
    data.get_temperature(), 
    data.get_energy_density(), 
    label=r'Set A', 
    color=color, 
    linewidth=2, 
    linestyle=linestyle
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\epsilon\, [\mathrm{GeV}^4]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.500
ymin = 0.0
ymax = 0.6
x_num_ticks = 7
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'energy_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color, linestyle
del fig, ax, xmin, xmax, ymin, ymax

####################################################################################################
# set A
print("Building plot: energy density, euler equation versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# e = T*s − P + sum_i ​mu_i rho_i
# e = energy density
# T = temperature 
# s = entropy density 
# P = pressure
# mu_i = i-th chemical potential​
# rho_i = i-th particle density
energy_density_euler = data.get_temperature()*data.get_entropy_density() - data.get_pressure()

# Plot data
color_e = 'black'
linestyle_e = '-'

color_euler = 'red'
linestyle_euler = '--'

ax.plot(
    data.get_temperature(), 
    data.get_energy_density(), 
    label=r'$\epsilon$', 
    color=color_e, 
    linewidth=2, 
    linestyle=linestyle_e
)

ax.plot(
    data.get_temperature(), 
    energy_density_euler, 
    label=r'$\epsilon = -P + T \, s + \sum_i \, \mu_i \rho_i$', 
    color=color_euler, 
    linewidth=2, 
    linestyle=linestyle_euler
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\epsilon\, [\mathrm{GeV}^4]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.500
ymin = 0.0
ymax = 0.6
x_num_ticks = 7
y_num_ticks = 5
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.68; auxY = 0.05
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'energy_euler_eq_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color_e, linestyle_e
del color_euler, linestyle_euler
del fig, ax, xmin, xmax, ymin, ymax
