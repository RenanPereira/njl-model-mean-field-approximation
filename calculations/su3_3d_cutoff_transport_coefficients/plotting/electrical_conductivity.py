import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.electrical_conductivity_data import ElectricalConductivityData


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_sigmae_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    ratio_data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    x_num_ticks: int = 6,
    y_num_ticks: int = 6,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.1f",
    annotation_texts: list[str] | None = None,
    x_annotation: float = 0.05,
    y_annotation: float = 0.05,
) -> tuple[Figure, Axes]:
    """
    ratio_data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_sigmae, label, color, linewidth, linestyle)
    """
    print("Building plot: electrical conductivity versus temperature.")

    print("Using datafiles:")
    for path_file_sigmae, _, _, _, _ in ratio_data_specs:
        print(path_file_sigmae)
    print()
    
    # Verify that the data provided have the same temperature grid
    datasets = []
    for path_file_sigmae, label, color, linewidth, linestyle in ratio_data_specs:
        data_sigmae = ElectricalConductivityData(path_file_sigmae)
        datasets.append((data_sigmae, label, color, linewidth, linestyle))

    for data_sigmae, label, color, linewidth, linestyle in datasets :
        if not np.array_equal(datasets[0][0].get_temperature(), data_sigmae.get_temperature()):
            raise ValueError("Temperature grids between datasets do not match.")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_sigmae, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data_sigmae.get_temperature(), 
            data_sigmae.get_electrical_conductivity(), 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )

    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=16, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\sigma_{\mathrm{e}}\, [\mathrm{GeV}]$', fontsize=20)

    # Configure axes using the helper function
    xmin = xlim[0]
    xmax = xlim[1]
    ymin = ylim[0]
    ymax = ylim[1]
    configure_axes(
        ax, 
        xmin, 
        xmax, 
        ymin, 
        ymax, 
        x_num_ticks, 
        y_num_ticks, 
        tick_fontsize=16, 
        spine_width=1.5, 
        tick_width=1.5, 
        tick_length=6
    )
    ax.xaxis.set_major_formatter(FormatStrFormatter(x_formatter))
    ax.yaxis.set_major_formatter(FormatStrFormatter(y_formatter))

    # Add text annotations
    if annotation_texts is not None:
        add_annotation_block(
            ax, 
            xmin, 
            xmax, 
            ymin, 
            ymax, 
            x_annotation, 
            y_annotation, 
            auxH=0.06, 
            texts=annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)

    return fig, ax


def plot_sigmae_over_temp_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    ratio_data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    x_num_ticks: int = 6,
    y_num_ticks: int = 6,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.1f",
    annotation_texts: list[str] | None = None,
    x_annotation: float = 0.05,
    y_annotation: float = 0.05,
) -> tuple[Figure, Axes]:
    """
    ratio_data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_sigmae, label, color, linewidth, linestyle)
    """
    print("Building plot: electrical conductivity over temperature versus temperature.")

    print("Using datafiles:")
    for path_file_sigmae, _, _, _, _ in ratio_data_specs:
        print(path_file_sigmae)
    print()
    
    # Verify that the data provided have the same temperature grid
    datasets = []
    for path_file_sigmae, label, color, linewidth, linestyle in ratio_data_specs:
        data_sigmae = ElectricalConductivityData(path_file_sigmae)
        datasets.append((data_sigmae, label, color, linewidth, linestyle))

    for data_sigmae, label, color, linewidth, linestyle in datasets :
        if not np.array_equal(datasets[0][0].get_temperature(), data_sigmae.get_temperature()):
            raise ValueError("Temperature grids between datasets do not match.")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_sigmae, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data_sigmae.get_temperature(), 
            data_sigmae.get_electrical_conductivity()/data_sigmae.get_temperature(), 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )

    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=16, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\sigma_{\mathrm{e}} / T$', fontsize=20)

    # Configure axes using the helper function
    xmin = xlim[0]
    xmax = xlim[1]
    ymin = ylim[0]
    ymax = ylim[1]
    configure_axes(
        ax, 
        xmin, 
        xmax, 
        ymin, 
        ymax, 
        x_num_ticks, 
        y_num_ticks, 
        tick_fontsize=16, 
        spine_width=1.5, 
        tick_width=1.5, 
        tick_length=6
    )
    ax.xaxis.set_major_formatter(FormatStrFormatter(x_formatter))
    ax.yaxis.set_major_formatter(FormatStrFormatter(y_formatter))

    # Add text annotations
    if annotation_texts is not None:
        add_annotation_block(
            ax, 
            xmin, 
            xmax, 
            ymin, 
            ymax, 
            x_annotation, 
            y_annotation, 
            auxH=0.06, 
            texts=annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)

    return fig, ax


##########################################################################


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"

ratio_datasets = [
    (
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV.dat",  
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY.dat",  
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_sigmae_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    ratio_datasets,
    path_output_plot_folder + "sigmae_vs_temp_CP0_setA.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.03),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.2f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu = 0.0\ \mathrm{GeV}$",
    ],
    x_annotation=0.03,
    y_annotation=0.6,
)

plot_sigmae_over_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    ratio_datasets,
    path_output_plot_folder + "sigmae_over_temp_vs_temp_CP0_setA.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.00, 0.15),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.2f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu = 0.0\ \mathrm{GeV}$",
    ],
    x_annotation=0.03,
    y_annotation=0.05,
)
