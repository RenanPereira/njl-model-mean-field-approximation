import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.io_utils import print_unique_filepaths
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData


# Select font that will be used for the different plots
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False
})


def plot_quark_masses_vs_temperature(    
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_effective_quark_mass: str | None = None,
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
    This function plots effective quark masses vs temperature from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: effective quark masses versus temperature.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, str, int, str]] = []
    for filepath, quark_species, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, quark_species, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, quark_species, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_temperature(),
            data.get_quark_effective_mass(quark_species), 
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
    if label_effective_quark_mass is None:
        ax.set_ylabel(r'$M_q\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_effective_quark_mass, fontsize=labels_fontsize)

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

def plot_normalized_quark_masses_vs_temperature(    
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_normalized_effective_quark_mass: str | None = None,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    x_num_ticks: int = 6,
    y_num_ticks: int = 6,
    tick_fontsize: int = 20,
    y_tick_max: float = 1.0,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.1f",
    annotation_texts: list[str] | None = None,
    x_annotation: float = 0.05,
    y_annotation: float = 0.05,
    annotation_vert_space: float = 0.06,
    annotation_fontsize: int = 18
) -> tuple[Figure, Axes]:
    """
    This function plots effective quark masses vs temperature from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: effective quark masses versus temperature.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, str, int, str]] = []
    for filepath, quark_species, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, quark_species, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, quark_species, label, color, linewidth, linestyle in datasets:
        quark_mass_vac = data.get_quark_effective_mass(quark_species)[0]
        
        ax.plot(
            data.get_temperature(),
            data.get_quark_effective_mass(quark_species)/quark_mass_vac, 
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
    if label_normalized_effective_quark_mass is None:
        ax.set_ylabel(r'$M_\ell/M_\ell^{\mathrm{vac}}$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_normalized_effective_quark_mass, fontsize=labels_fontsize)

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
        tick_length=6,
        y_tick_max=y_tick_max
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

def plot_quark_masses_vs_baryon_chem_pot(    
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_effective_quark_mass: str | None = None,
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
    This function plots effective quark masses vs baryon chemical potential from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: effective quark masses versus baryon chemical potential.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, str, int, str]] = []
    for filepath, quark_species, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, quark_species, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, quark_species, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_baryon_chemical_potential(),
            data.get_quark_effective_mass(quark_species), 
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
    ax.set_xlabel(r'$\mu_B\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    if label_effective_quark_mass is None:
        ax.set_ylabel(r'$M_q\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_effective_quark_mass, fontsize=labels_fontsize)

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
