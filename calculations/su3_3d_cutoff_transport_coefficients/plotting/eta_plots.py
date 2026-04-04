import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy.interpolate import interp1d
import numpy as np
import math
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.shear_viscosity_data import ShearViscosityData
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_eta_vs_temp(
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
) -> tuple[Figure, Axes]:
    """
    data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_eta, label, color, linewidth, linestyle)
    """
    print("Building plot: shear viscosity versus temperature.")
    
    print("Using datafiles:")
    for path_file_eta, _, _, _, _ in data_specs:
        print(path_file_eta)
    print()

    datasets: list[tuple[ShearViscosityData, str, str, int, str]] = []
    for path_file_eta, label, color, linewidth, linestyle in data_specs:
        data_eta = ShearViscosityData(path_file_eta)
        datasets.append((data_eta, label, color, linewidth, linestyle))
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_eta, label, color, linewidth, linestyle in datasets:
        temp = data_eta.get_temperature()
        eta = data_eta.get_shear_viscosity()
        ax.plot(
            temp, 
            eta, 
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
    ax.set_ylabel(r'$\eta\, [\mathrm{GeV}^3]$', fontsize=20)

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


def plot_eta_over_s_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_file_thermodynamics: str,
    path_output_plot: str,
    include_kss_bound: bool = True,
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
        (path_file_eta, label, color, linewidth, linestyle)
    """
    print("Building plot: shear viscosity over entropy density ratio versus temperature.")
    
    print("Using datafiles:")
    for path_file_eta, _, _, _, _ in data_specs:
        print(path_file_eta)
    print()

    datasets: list[tuple[ShearViscosityData, str, str, int, str]] = []
    for path_file_eta, label, color, linewidth, linestyle in data_specs:
        data_eta = ShearViscosityData(path_file_eta)
        datasets.append((data_eta, label, color, linewidth, linestyle))

    # Get entropy from thermodynamics data (for this parameter set) and interpolate it
    data_thermodynamics = FixedChemPotTempData(path_file_thermodynamics)
    entropy_dens_interpolation = interp1d(
        data_thermodynamics.get_temperature(), 
        data_thermodynamics.get_entropy_density(), 
        kind='linear'
    )
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data_eta, label, color, linewidth, linestyle in datasets:
        temp = data_eta.get_temperature()
        eta = data_eta.get_shear_viscosity()
        entropy_dens = entropy_dens_interpolation(temp)
        ax.plot(
            temp, 
            eta/entropy_dens, 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )

    if include_kss_bound:
        # Kovtun-Son-Starinets conformal limit (KSS), eta/s|KSS = 1/4*pi
        ax.axhline(
            y=1.0/(4.0*math.pi),
            color='black',
            linestyle='--',
            linewidth=2,
            label=r'KSS'
        )
    
    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=16, frameon=False)
        
    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\eta/s$', fontsize=20)

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
