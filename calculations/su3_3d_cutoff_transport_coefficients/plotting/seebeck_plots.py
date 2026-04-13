import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np

from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.seebeck_sigmae_product_data import SeebeckSigmaeProductData
from common_utils.electrical_conductivity_data import ElectricalConductivityData


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_seebeck_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
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
    data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_seebeck_sigmae_prod, path_file_sigmae, label, color, linewidth, linestyle)
    """
    print("Building plot: seebeck versus temperature.")   
    
    print("Using datafiles:")
    for path_file_seebeck_sigmae_prod, path_file_sigmae, _, _, _, _ in data_specs:
        print(path_file_seebeck_sigmae_prod)
        print(path_file_sigmae)
    print()

    datasets: list[tuple[SeebeckSigmaeProductData, ElectricalConductivityData, str, str, int, str]] = []
    for path_file_seebeck_sigmae_prod, path_file_sigmae, label, color, linewidth, linestyle in data_specs:
        data_seebeck_sigmae_prod = SeebeckSigmaeProductData(path_file_seebeck_sigmae_prod)
        data_sigmae = ElectricalConductivityData(path_file_sigmae)
        datasets.append((data_seebeck_sigmae_prod, data_sigmae, label, color, linewidth, linestyle))

    # Verify that the seebeck_sigmae_prod and sigmae paired data provided have the same temperature grid
    for data_seebeck_sigmae_prod, data_sigmae, label, color, linewidth, linestyle in datasets :
        if not np.array_equal(data_seebeck_sigmae_prod.get_temperature(), data_sigmae.get_temperature()):
            raise ValueError("Temperature grids between paired seebeck_sigmae_prod and sigmae do not match.")
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_seebeck_sigmae_prod, data_sigmae, label, color, linewidth, linestyle in datasets:
        temp = data_seebeck_sigmae_prod.get_temperature()
        seebeck_sigmae = data_seebeck_sigmae_prod.get_seebeck_sigmae_product()
        sigma_e = data_sigmae.get_electrical_conductivity()
        ax.plot(
            temp, 
            seebeck_sigmae/sigma_e, 
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
    ax.set_xlabel(r"$T\, [\mathrm{GeV}]$", fontsize=20)
    ax.set_ylabel(r"$\mathrm{S}$", fontsize=20)

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
