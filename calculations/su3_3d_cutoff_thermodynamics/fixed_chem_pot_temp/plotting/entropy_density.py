import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from common_utils.io_utils import print_unique_filepaths
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from common_utils.stefan_boltzmann import StefanBoltzmannMasslessFlavorDegenerateQuarks


# Select font that will be used for the different plots
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False
})


def plot_entropy_density_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_entropy_density: str | None = None,
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
    This function plots the entropy density vs temperature from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, label, color, linewidth, linestyle)
    """
    print("Building plot: effective entropy density versus temperature.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str]] = []
    for filepath, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, label, color, linewidth, linestyle in datasets:
        ax.plot(
            data.get_temperature(),
            data.get_entropy_density(),
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
    if label_entropy_density is None:
        ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_entropy_density, fontsize=labels_fontsize)

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

def plot_entropy_density_dPdT_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str, bool]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_entropy_density: str | None = None,
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
    This function plots the entropy density, dPdT vs temperature from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, label, color, linewidth, linestyle)
    """
    print("Building plot: effective entropy density, dPdT versus temperature.")
    print_unique_filepaths(data_specs)

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str, bool]] = []
    for filepath, label, color, linewidth, linestyle, dPdT_calculation in data_specs:
        data = FixedChemPotTempData(filepath)
        datasets.append((data, label, color, linewidth, linestyle, dPdT_calculation))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, label, color, linewidth, linestyle, dPdT_calculation in datasets:
        if (not dPdT_calculation):
            ax.plot(
                data.get_temperature(),
                data.get_entropy_density(),  
                label=label, 
                color=color, 
                linewidth=linewidth, 
                linestyle=linestyle
            )
        else:
            dPdT = np.gradient(data.get_pressure(), data.get_temperature())
            ax.plot(
                data.get_temperature(),
                dPdT,  
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
    if label_entropy_density is None:
        ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_entropy_density, fontsize=labels_fontsize)

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

def plot_s_over_temp3_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    stefan_boltzmann_limit: bool = True,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
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
    data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_eta, label, color, linewidth, linestyle)
    """
    print("Building plot: entropy density over temperature^3 versus temperature.")

    print("Using datafiles:")
    for path_file_s, _, _, _, _ in data_specs:
        print(path_file_s)
    print()

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str]] = []
    for path_file_s, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(path_file_s)
        datasets.append((data, label, color, linewidth, linestyle))

    # Verify that the data provided have the same temperature grid
    for data, label, color, linewidth, linestyle in datasets :
        if not np.array_equal(datasets[0][0].get_temperature(), data.get_temperature()):
            raise ValueError("Temperature grids between datasets do not match.")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, label, color, linewidth, linestyle in datasets:
        
        # remove zero temperature values
        temp = data.get_temperature()
        s = data.get_entropy_density()
        mask = temp > 0
        temp = temp[mask]
        s = s[mask]
        """
        temp = []
        s = []
        for i in range(0, data.size()):
            if data.get_temperature()[i]>0:
                temp.append(data.get_temperature()[i])
                s.append(data.get_entropy_density()[i]) 
        temp = np.array(temp)
        s = np.array(s)
        """

        ax.plot(
            temp, 
            s/(temp**3), 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )
    
    # add Stefan Boltzmann reference
    if stefan_boltzmann_limit:
        sb = StefanBoltzmannMasslessFlavorDegenerateQuarks(number_of_colors=3,number_of_flavors=3)
        temp = datasets[0][0].get_temperature()
        mask = temp > 0
        temp = temp[mask]
        ax.plot(
            temp, 
            sb.entropy_density_over_temp3(0.0, temp), 
            label=r"SB", 
            color="Black", 
            linewidth=2, 
            linestyle="--"
        )

    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    ax.set_ylabel(r'$s / T^3$', fontsize=labels_fontsize)
    
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
            auxH=annotation_vert_space, 
            texts=annotation_texts, 
            fontsize=annotation_fontsize
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)

    return fig, ax

def plot_entropy_density_vs_baryon_chem_pot(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
    label_entropy_density: str | None = None,
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
    This function plots entropy_density vs baryon chemical potential from FixedChemPotTempData.
    
    data_specs:
        List of tuples defining datasets and plot styles:
        (filepath, quark_species, label, color, linewidth, linestyle)
    """
    print("Building plot: entropy_density versus baryon chemical potential.")
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
            data.get_entropy_density(),
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
    if label_entropy_density is None:
        ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=labels_fontsize)
    else: 
        ax.set_ylabel(label_entropy_density, fontsize=labels_fontsize)

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

def plot_s_over_temp3_vs_baryon_chem_pot(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    stefan_boltzmann_limit: bool = True,
    legend_loc: str | None = None,
    legend_fontsize: int = 18,
    labels_fontsize: int = 22,
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
    data_specs:
        List of tuples defining datasets and plot styles:
        (path_file_eta, label, color, linewidth, linestyle)
    """
    print("Building plot: entropy density over temperature^3 versus baryon chemical potential.")

    print("Using datafiles:")
    for path_file_s, _, _, _, _ in data_specs:
        print(path_file_s)
    print()

    datasets: list[tuple[FixedChemPotTempData, str, str, int, str]] = []
    for path_file_s, label, color, linewidth, linestyle in data_specs:
        data = FixedChemPotTempData(path_file_s)
        datasets.append((data, label, color, linewidth, linestyle))

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    for data, label, color, linewidth, linestyle in datasets:
        
        # remove zero temperature values
        temp = data.get_temperature()
        s = data.get_entropy_density()
        mask = temp > 0
        temp = temp[mask]
        s = s[mask]
        """
        temp = []
        s = []
        for i in range(0, data.size()):
            if data.get_temperature()[i]>0:
                temp.append(data.get_temperature()[i])
                s.append(data.get_entropy_density()[i]) 
        temp = np.array(temp)
        s = np.array(s)
        """

        ax.plot(
            data.get_baryon_chemical_potential(), 
            s/(temp**3), 
            label=label, 
            color=color, 
            linewidth=linewidth, 
            linestyle=linestyle
        )
    
    # add Stefan Boltzmann reference using the first dataset provided
    if stefan_boltzmann_limit:
        sb = StefanBoltzmannMasslessFlavorDegenerateQuarks(number_of_colors=3,number_of_flavors=3)
        temp = datasets[0][0].get_temperature()
        mask = temp > 0
        temp = temp[mask]
        cpu = datasets[0][0].get_quark_chemical_potential("up_quark")
        ax.plot(
            cpu, 
            sb.entropy_density_over_temp3(cpu, temp), 
            label=r"SB", 
            color="Black", 
            linewidth=2, 
            linestyle="--"
        )

    # Grid
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Legend
    if legend_loc is not None:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$\mu_B\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    ax.set_ylabel(r'$s / T^3$', fontsize=labels_fontsize)
    
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
            auxH=annotation_vert_space, 
            texts=annotation_texts, 
            fontsize=annotation_fontsize
        )

    fig.tight_layout()

    plt.savefig(path_output_plot)

    return fig, ax
