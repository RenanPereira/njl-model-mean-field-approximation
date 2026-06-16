import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from common_utils.io_utils import print_unique_filepaths
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData

# Select font that will be used for the different plots
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False
})


def plot_pressure_vs_energy(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_energy_density: str | None = None,
    label_pressure: str | None = None,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    x_num_ticks: int = 6,
    y_num_ticks: int = 6,
    tick_fontsize: int = 20,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.1f",
    annotation_texts: list[str] | None = None,
    x_annotation: float = 0.05,
    y_annotation: float = 0.05,
    annotation_vert_space: float = 0.06,
    annotation_fontsize: int = 18
) -> tuple[Figure, Axes]:
    """
    This function plots the pressure versus energy density from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, label, color, linewidth, linestyle)
    """
    print("Building plot: pressure versus energy density.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str]] = []
    for filepath, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_energy_density(),
            data.get_pressure(),  
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )
        
    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize, frameon=False)

    # Axes labels
    if label_energy_density is None:
        ax.set_xlabel(r'$\epsilon \, [\mathrm{GeV}^4]$', fontsize=labels_fontsize)
    else:
        ax.set_xlabel(label_energy_density, fontsize=labels_fontsize)
    if label_pressure is None:
        ax.set_ylabel(r'$P \, [\mathrm{GeV}^4]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_pressure, fontsize=labels_fontsize)

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
        tick_fontsize=tick_fontsize, 
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
            annotation_vert_space, 
            annotation_texts, 
            fontsize=annotation_fontsize
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)
    
    return fig, ax
