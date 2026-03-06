import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.electrical_conductivity_data import ElectricalConductivityData


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_plots_folder = "su3_3d_cutoff_transport_coefficients/plots/"


####################################################################################################
# set A
print("Building plot: electrical conductivity as a function of temperature.")
print("Parameter set A, zero chemical potential and different methods for the integrated cross section.\n")

parameter_set = "setA"

method = "COMPLETE_COV"
data_complete = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

method = "KLEVANSKY"
data_klevansky = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

method = "ZHUANG"
data_zhuang = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_complete = 'black'
linestyle_complete = '-'

color_klevansky = 'red'
linestyle_klevansky = '-'

color_zhuang = 'blue'
linestyle_zhuang = '-'


ax.plot(
    data_complete.get_temperature(), 
    data_complete.get_electrical_conductivity(), 
    label=r'Method I', 
    color=color_complete, 
    linewidth=2, 
    linestyle=linestyle_complete
)

ax.plot(
    data_klevansky.get_temperature(), 
    data_klevansky.get_electrical_conductivity(), 
    label=r'Method II', 
    color=color_klevansky, 
    linewidth=2, 
    linestyle=linestyle_klevansky
)

ax.plot(
    data_zhuang.get_temperature(), 
    data_zhuang.get_electrical_conductivity(), 
    label=r'Method III', 
    color=color_zhuang, 
    linewidth=2, 
    linestyle=linestyle_zhuang
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\mathrm{e}}\, [\mathrm{GeV}]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=14, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.120
xmax = 0.300
ymin = 0.00
ymax = 0.03
x_num_ticks = 4
y_num_ticks = 4
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.03; auxY = 0.6
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Automatically adjust layout

fig.tight_layout()

plotname = f'electrical_conductivity_CP0_{parameter_set}.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, method, plotname
del data_complete, data_klevansky, data_zhuang
del color_complete, linestyle_complete
del color_klevansky, linestyle_klevansky
del color_zhuang, linestyle_zhuang
del fig, ax, xmin, xmax, ymin, ymax


####################################################################################################
# set A
print("Building plot: electrical conductivity to temperature ratio as a function of temperature.")
print("Parameter set A, zero chemical potential and different methods for the integrated cross section.\n")

parameter_set = "setA"

method = "COMPLETE_COV"
data_complete = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

method = "KLEVANSKY"
data_klevansky = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

method = "ZHUANG"
data_zhuang = ElectricalConductivityData(path_data_folder + f'ElectricalConductivity_{parameter_set}_{method}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_complete = 'black'
linestyle_complete = '-'

color_klevansky = 'red'
linestyle_klevansky = '-'

color_zhuang = 'blue'
linestyle_zhuang = '-'


ax.plot(
    data_complete.get_temperature(), 
    data_complete.get_electrical_conductivity()/data_complete.get_temperature(), 
    label=r'Method I', 
    color=color_complete, 
    linewidth=2, 
    linestyle=linestyle_complete
)

ax.plot(
    data_klevansky.get_temperature(), 
    data_klevansky.get_electrical_conductivity()/data_klevansky.get_temperature(), 
    label=r'Method II', 
    color=color_klevansky, 
    linewidth=2, 
    linestyle=linestyle_klevansky
)

ax.plot(
    data_zhuang.get_temperature(), 
    data_zhuang.get_electrical_conductivity()/data_zhuang.get_temperature(), 
    label=r'Method III', 
    color=color_zhuang, 
    linewidth=2, 
    linestyle=linestyle_zhuang
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\sigma_{\mathrm{e}} / T$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper left", fontsize=14, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.120
xmax = 0.300
ymin = 0.00
ymax = 0.15
x_num_ticks = 4
y_num_ticks = 6
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add text annotations
auxH = 0.06; auxX = 0.03; auxY = 0.05
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Automatically adjust layout

fig.tight_layout()

plotname = f'electrical_conductivity_temp_ratio_CP0_{parameter_set}.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, method, plotname
del data_complete, data_klevansky, data_zhuang
del color_complete, linestyle_complete
del color_klevansky, linestyle_klevansky
del color_zhuang, linestyle_zhuang
del fig, ax, xmin, xmax, ymin, ymax
