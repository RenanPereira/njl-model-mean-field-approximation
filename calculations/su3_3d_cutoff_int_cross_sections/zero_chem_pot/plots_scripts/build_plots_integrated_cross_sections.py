import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import configure_axes
from common_utils.integrated_cross_section_data import IntegratedCrossSectionData


PROCESS_YLABELS = {
    "UUUU": r'$\overline{\sigma}_{uu \rightarrow uu} \, [\mathrm{GeV}^{-2}]$',
    "UUBarUUBar": r'$\overline{\sigma}_{u\bar{u} \rightarrow u\bar{u}} \, [\mathrm{GeV}^{-2}]$',
    "UUBarDDBar": r'$\overline{\sigma}_{u\bar{u} \rightarrow d\bar{d}} \, [\mathrm{GeV}^{-2}]$',
    "UUBarSSBar": r'$\overline{\sigma}_{u\bar{u} \rightarrow s\bar{s}} \, [\mathrm{GeV}^{-2}]$',
    "UDUD": r'$\overline{\sigma}_{ud \rightarrow ud} \, [\mathrm{GeV}^{-2}]$',
    "USUS": r'$\overline{\sigma}_{us \rightarrow us} \, [\mathrm{GeV}^{-2}]$',
    "USBarUSBar": r'$\overline{\sigma}_{u\bar{s} \rightarrow u\bar{s}} \, [\mathrm{GeV}^{-2}]$',
    "UDBarUDBar": r'$\overline{\sigma}_{u\bar{d} \rightarrow u\bar{d}} \, [\mathrm{GeV}^{-2}]$',
    "SSSS": r'$\overline{\sigma}_{ss \rightarrow ss} \, [\mathrm{GeV}^{-2}]$',
    "SSBarUUBar": r'$\overline{\sigma}_{s\bar{s} \rightarrow u\bar{u}} \, [\mathrm{GeV}^{-2}]$',
    "SSBarSSBar": r'$\overline{\sigma}_{s\bar{s} \rightarrow s\bar{s}} \, [\mathrm{GeV}^{-2}]$',
}


def process_to_ylabel_latex(process: str) -> str:
    try:
        return PROCESS_YLABELS[process]
    except KeyError:
        raise ValueError(f"Unknown process '{process}'")


def plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    set,
    process,
    legend_loc,
    xmin,
    xmax,
    ymin,
    ymax,
    x_num_ticks,
    y_num_ticks
):
    complete_label = "Complete"
    klevansky_label = "Klevansky"
    zhuang_label = "Zhuang"

    print(f'Building plot: integrated cross section as a function of temperature for {process} (zero chemical potential)')

    # Load data from the file
    filename_complete  = f'IntegratedCrossSection_{set}_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat'
    data_complete = IntegratedCrossSectionData(data_folder + filename_complete )

    filename_klevansky  = f'IntegratedCrossSection_{set}_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat'
    data_klevansky = IntegratedCrossSectionData(data_folder + filename_klevansky )

    filename_zhuang  = f'IntegratedCrossSection_{set}_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat'
    data_zhuang = IntegratedCrossSectionData(data_folder + filename_zhuang )

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    # Plot data
    color_complete = 'black'
    linestyle_complete = '-'
    color_klevansky = 'red'
    linestyle_klevansky = '-'
    color_zhuang = 'blue'
    linestyle_zhuang = '-'

    # Plot cross sections
    ax.plot(data_complete.get_temperature(), data_complete.get_integrated_cross_section(), label=complete_label, color=color_complete, linewidth=2, linestyle=linestyle_complete)
    ax.plot(data_klevansky.get_temperature(), data_klevansky.get_integrated_cross_section(), label=klevansky_label, color=color_klevansky, linewidth=2, linestyle=linestyle_klevansky)
    ax.plot(data_zhuang.get_temperature(), data_zhuang.get_integrated_cross_section(), label=zhuang_label, color=color_zhuang, linewidth=2, linestyle=linestyle_zhuang)

    # # Axes labels
    ax.set_xlabel(r'T$\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(process_to_ylabel_latex(process), fontsize=20)

    # # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc=legend_loc, fontsize=14, frameon=False, title='Method', title_fontsize=14)

    # # Configure axes using the helper function
    configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # # Automatically adjust layout
    fig.tight_layout()

    plotname = f'integrated_cross_section_{set}_{process.lower()}_CP0.png'
    plt.savefig(plots_folder + plotname)

    # Clean up
    plt.clf()
    plt.close()


#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

# Common configurations between plots
fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
data_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/data/"
plots_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/plots/"


plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UUUU",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UUBarUUBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    12.0,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UUBarDDBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UUBarSSBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    3,
    4,
    4
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "USUS",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "USBarUSBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UDUD",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "UDBarUDBar",
    'upper left',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "SSSS",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "SSBarUUBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    8,
    4,
    5
)

plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_folder,
    plots_folder,
    "setA",
    "SSBarSSBar",
    'upper right',
    0.120,
    0.300,
    0.0,
    12.0,
    4,
    5
)
