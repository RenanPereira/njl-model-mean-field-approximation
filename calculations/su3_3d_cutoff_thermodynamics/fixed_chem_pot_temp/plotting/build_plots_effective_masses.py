import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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
print("Building plot: effective masses versus temperature.")
print("Parameter set A, zero chemical potential.\n")

parameter_set = "setA"
suffix = "TMin0p000000_TMax0p500000_CP0"
data = FixedChemPotTempData(path_data_folder + f'SU3NJL3DCutoffFixedChemPotTemp_{parameter_set}_{suffix}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_l = 'black'
linestyle_l = '-'

color_s = 'red'
linestyle_s = '-'

ax.plot(
    data.get_temperature(), 
    data.get_up_quark_effective_mass(), 
    label=r'$M_l$', 
    color=color_l, 
    linewidth=2, 
    linestyle=linestyle_l
)

ax.plot(
    data.get_temperature(), 
    data.get_strange_quark_effective_mass(), 
    label=r'$M_s$', 
    color=color_s, 
    linewidth=2, 
    linestyle=linestyle_s
)

# # Axes labels
ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$M_q\, [\mathrm{GeV}]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper right", fontsize=16, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.000
xmax = 0.500
ymin = 0.0
ymax = 0.6
x_num_ticks = 6
y_num_ticks = 7
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add text annotations
auxH = 0.06; auxX = 0.05; auxY = 0.05
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

fig.tight_layout()

plotname = f'quark_eff_masses_vs_temp_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, plotname
del data
del color_l, linestyle_l, color_s, linestyle_s
del fig, ax, xmin, xmax, ymin, ymax
