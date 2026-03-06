import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import configure_axes
from common_utils.integrated_cross_section_data import IntegratedCrossSectionData


#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


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
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_folder: str,
    plots_folder: str,
    parameter_set: str,
    process: str,
    legend_loc: str,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    x_num_ticks: int,
    y_num_ticks: int
) -> None:
    complete_label = "Method I"
    klevansky_label = "Method II"
    zhuang_label = "Method III"

    print(f'Building plot: integrated cross section as a function of temperature for {process} (zero chemical potential)')

    # Load data from the file
    data_complete = IntegratedCrossSectionData.from_matching_files(data_folder, parameter_set, process, "COMPLETE_COV")
    data_klevansky = IntegratedCrossSectionData.from_matching_files(data_folder, parameter_set, process, "KLEVANSKY")
    data_zhuang = IntegratedCrossSectionData.from_matching_files(data_folder, parameter_set, process, "ZHUANG")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    # Plot data
    color_complete = 'black'
    linestyle_complete = '-'
    color_klevansky = 'red'
    linestyle_klevansky = '-'
    color_zhuang = 'blue'
    linestyle_zhuang = '-'

    # Plot integrated cross sections
    ax.plot(data_complete.get_temperature(), data_complete.get_integrated_cross_section(), label=complete_label, color=color_complete, linewidth=2, linestyle=linestyle_complete)
    ax.plot(data_klevansky.get_temperature(), data_klevansky.get_integrated_cross_section(), label=klevansky_label, color=color_klevansky, linewidth=2, linestyle=linestyle_klevansky)
    ax.plot(data_zhuang.get_temperature(), data_zhuang.get_integrated_cross_section(), label=zhuang_label, color=color_zhuang, linewidth=2, linestyle=linestyle_zhuang)

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(process_to_ylabel_latex(process), fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc=legend_loc, fontsize=14, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Automatically adjust layout
    fig.tight_layout()

    plotname = f'integrated_cross_section_{parameter_set}_{process.lower()}_CP0.png'
    plt.savefig(plots_folder + plotname)

    # Clean up
    plt.clf()
    plt.close()


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
    "upper right",
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
    "upper right",
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
    "upper right",
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
    "upper right",
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
    "upper right",
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
    "upper right",
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
    "upper right",
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
    "upper left",
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
    "upper right",
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
    "upper right",
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
    "upper right",
    0.120,
    0.300,
    0.0,
    12.0,
    4,
    5
)
