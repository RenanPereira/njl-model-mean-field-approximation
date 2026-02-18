import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import *
from common_utils.quark_relaxation_times_data import QuarkRelaxationTimesData


####################################################################################################
# Common configurations between plots

#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_plots_folder = "su3_3d_cutoff_quark_relaxation_times/plots/"


####################################################################################################
# set A
print("Building plot: quark relaxation time with set A, at zero chemical potential and using different methods for the integrated cross section")

parameter_set = "setA"

method = "COMPLETE_COV"
data_complete = QuarkRelaxationTimesData(path_data_folder + f'RelaxationTimes_{parameter_set}_{method}.dat')

method = "KLEVANSKY"
data_klevansky = QuarkRelaxationTimesData(path_data_folder + f'RelaxationTimes_{parameter_set}_{method}.dat')

method = "ZHUANG"
data_zhuang = QuarkRelaxationTimesData(path_data_folder + f'RelaxationTimes_{parameter_set}_{method}.dat')

# Create a new figure
fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

# Plot data
color_complete_u = 'black'
linestyle_complete_u = '-'
color_complete_s = 'red'
linestyle_complete_s = '-'

color_klevansky_u = 'black'
linestyle_klevansky_u = '--'
color_klevansky_s = 'red'
linestyle_klevansky_s = '--'

color_zhuang_u = 'black'
linestyle_zhuang_u = ':'
color_zhuang_s = 'red'
linestyle_zhuang_s = ':'

# Plot relaxation times
ax.plot(data_complete.get_temperature(), data_complete.get_rel_time_up_quark(), label=r'$\tau_{l,\overline{l}}$  | Method I', color=color_complete_u, linewidth=2, linestyle=linestyle_complete_u)
ax.plot(data_complete.get_temperature(), data_complete.get_rel_time_strange_quark(), label=r'$\tau_{s,\overline{s}}$ | Method I', color=color_complete_s, linewidth=2, linestyle=linestyle_complete_s)

ax.plot(data_klevansky.get_temperature(), data_klevansky.get_rel_time_up_quark(), label=r'$\tau_{l,\overline{l}}$  | Method II', color=color_klevansky_u, linewidth=2, linestyle=linestyle_klevansky_u)
ax.plot(data_klevansky.get_temperature(), data_klevansky.get_rel_time_strange_quark(), label=r'$\tau_{s,\overline{s}}$ | Method II', color=color_klevansky_s, linewidth=2, linestyle=linestyle_klevansky_s)

ax.plot(data_zhuang.get_temperature(), data_zhuang.get_rel_time_up_quark(), label=r'$\tau_{l,\overline{l}}$  | Method III', color=color_zhuang_u, linewidth=2, linestyle=linestyle_zhuang_u)
ax.plot(data_zhuang.get_temperature(), data_zhuang.get_rel_time_strange_quark(), label=r'$\tau_{s,\overline{s}}$ | Method III', color=color_zhuang_s, linewidth=2, linestyle=linestyle_zhuang_s)

# # Axes labels
ax.set_xlabel(r'T$\, [\mathrm{GeV}]$', fontsize=20)
ax.set_ylabel(r'$\tau\, [\mathrm{fm}]$', fontsize=20)

# # Grid and legend
ax.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc="upper right", fontsize=14, frameon=False, title_fontsize=14)

# # Configure axes using the helper function
xmin = 0.120
xmax = 0.300
ymin = 0.0
ymax = 10
x_num_ticks = 4
y_num_ticks = 6
configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Add text annotations
auxH = 0.06; auxX = 0.03; auxY = 0.03
texts = [
    r'set A',
    r'$\mu = 0.0\ \mathrm{GeV}$',
]
add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)


# Automatically adjust layout

fig.tight_layout()

plotname = f'quark_relaxation_time_{parameter_set}_CP0.png'
plt.savefig(path_plots_folder + plotname)

# Clean up
plt.clf()
plt.close()

# Clear variables (optional cleanup)
del parameter_set, method, plotname
del data_complete, data_klevansky, data_zhuang
del color_complete_u, linestyle_complete_u, color_complete_s, linestyle_complete_s
del color_klevansky_u, linestyle_klevansky_u, color_klevansky_s, linestyle_klevansky_s
del color_zhuang_u, linestyle_zhuang_u, color_zhuang_s, linestyle_zhuang_s
del fig, ax, xmin, xmax, ymin, ymax
