import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.integrated_cross_section_data import IntegratedCrossSectionData


#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_integrated_cross_section_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[list[str], str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    label_int_cross_section: str | None = None,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    x_num_ticks: int = 6,
    y_num_ticks: int = 6,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.1f",
    annotation_texts: list[str] | None = None,
    x_annotation: float = 0.05,
    y_annotation: float = 0.05,
    annotation_vert_space: float = 0.06
) -> tuple[Figure, Axes]:
    """
    This function plots integrated cross sections vs temperature.
    
    The variable `data_specs` provides all the necessary data and information for the plot.
    It is a tuple where each tuple entry is considered an independent dataset.
    
    The list of filepaths provided for each entry in the tuple must be related to the same integrated cross sections data.
    These data files are joined together when the class IntegratedCrossSectionData is inicialized for each entry in the tuple.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (list[filepath_int_cross_section], label, color, linewidth, linestyle)
    """
    print("Building plot: integrated cross section versus temperature")

    print("The following data files will be used (each set of datafiles is merged):")
    data_set_number = 1
    for list_filepath, _, _, _, _ in data_specs:
        print(f"Set {data_set_number} of datafiles:")
        for path in list_filepath:
            print(path)
        data_set_number += 1 
    print()
    
    datasets: list[tuple[IntegratedCrossSectionData, str, str, int, str]] = []
    
    for list_filepath, label, color, linewidth, linestyle in data_specs:
        data_int_cross_section = IntegratedCrossSectionData(list_filepath)
        datasets.append((data_int_cross_section, label, color, linewidth, linestyle))
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_int_cross_section, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data_int_cross_section.get_temperature(),
            data_int_cross_section.get_integrated_cross_section(), 
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
    if label_int_cross_section is None:
        ax.set_ylabel(r'$\overline{\sigma} \, [\mathrm{GeV}^{-2}]$', fontsize=20)
    else: 
        ax.set_ylabel(label_int_cross_section, fontsize=20)

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
            auxH=annotation_vert_space, 
            texts=annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)
    
    return fig, ax
