import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.quark_relaxation_times_data import QuarkRelaxationTimesData


# Select font that will be used for the different plots
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False
})


def plot_quark_rel_time_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_rel_time: str | None = None,
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
    This function plots quark relaxation time vs temperature.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath_relaxation_time, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: quark relaxation time versus temperature")
    
    datasets: list[tuple[QuarkRelaxationTimesData, str, str, str, int, str]] = []
    
    for filepath, quark_species, label, color, linewidth, linestyle in data_specs:
        data_quark_rel_time = QuarkRelaxationTimesData(filepath)
        datasets.append((data_quark_rel_time, quark_species, label, color, linewidth, linestyle))
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_quark_rel_time, quark_species, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data_quark_rel_time.get_temperature(),
            data_quark_rel_time.get_rel_time(quark_species), 
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
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    if label_rel_time is None:
        ax.set_ylabel(r'$\tau\, [\mathrm{fm}]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_rel_time, fontsize=labels_fontsize)
    
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

    fig.savefig(path_output_plot)
    
    return fig, ax
