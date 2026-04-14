import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from common_utils.io_utils import print_unique_filepaths
from common_utils.physical_constants import hbarc_gevfm


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_quark_dens_vs_chem_pot(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, str, int, str]],
    x_axis_chem_pot_quark_species: str,
    path_output_plot: str,
    legend_loc: str | None = None,
    label_quark_density: str | None = None,
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
    This function plots quark densities vs checmical potential from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: quark density versus chemical potential.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, str, int, str]] = []
    for filepath, quark_species, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, quark_species, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    for data, quark_species, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_quark_chemical_potential(x_axis_chem_pot_quark_species),
            data.get_quark_density(quark_species)/(hbarc_gevfm**3), 
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
    ax.set_xlabel(r'$\mu\, [\mathrm{GeV}]$', fontsize=20)
    if label_quark_density is None:
        ax.set_ylabel(r'$\rho_q\, [\mathrm{fm}^{-3}]$', fontsize=20)
    else: 
        ax.set_ylabel(label_quark_density, fontsize=20)

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
            annotation_vert_space, 
            annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)
    
    return fig, ax


def plot_baryon_dens_vs_baryon_chem_pot(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
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
    annotation_vert_space: float = 0.06
) -> tuple[Figure, Axes]:
    """
    This function plots the baryon density vs baryon checmical potential from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, label, color, linewidth, linestyle)
    """
    print("Building plot: quark density versus chemical potential.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str]] = []
    for filepath, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    for data, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_baryon_chemical_potential(),
            data.get_baryon_density()/(hbarc_gevfm**3), 
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
    ax.set_xlabel(r'$\mu_B\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\rho_B\, [\mathrm{fm}^{-3}]$', fontsize=20)

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
            annotation_vert_space, 
            annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)
    
    return fig, ax


def plot_baryon_dens_dPdmuB_vs_baryon_chem_pot(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str, str, str, int, str]],
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
    annotation_vert_space: float = 0.06
) -> tuple[Figure, Axes]:
    """
    This function plots the baryon density vs baryon checmical potential from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, label, color, linewidth, linestyle, label_dPdmuB, color_dPdmuB, linewidth_dPdmuB, linestyle_dPdmuB)
    """
    print("Building plot: baryon density and dP/dmuB versus baryon chemical potential.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str, str, str, int, str]] = []
    for filepath, label, color, linewidth, linestyle, label_dPdmuB, color_dPdmuB, linewidth_dPdmuB, linestyle_dPdmuB in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, label, color, linewidth, linestyle, label_dPdmuB, color_dPdmuB, linewidth_dPdmuB, linestyle_dPdmuB))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    for data, label, color, linewidth, linestyle, label_dPdmuB, color_dPdmuB, linewidth_dPdmuB, linestyle_dPdmuB in datasets:
        ax.plot(
            data.get_baryon_chemical_potential(),
            data.get_baryon_density()/(hbarc_gevfm**3), 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )

        dPdmuB = np.gradient(data.get_pressure(), data.get_baryon_chemical_potential())
        ax.plot(
            data.get_baryon_chemical_potential(),
            dPdmuB/(hbarc_gevfm**3), 
            label=label_dPdmuB, 
            color=color_dPdmuB, 
            linewidth=linewidth_dPdmuB, 
            linestyle=linestyle_dPdmuB
        )
    
    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=16, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$\mu_B\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\rho_B\, [\mathrm{fm}^{-3}]$', fontsize=20)

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
            annotation_vert_space, 
            annotation_texts, 
            fontsize=16
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)
    
    return fig, ax
