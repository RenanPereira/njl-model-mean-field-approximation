import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from common_utils.stefan_boltzmann import StefanBoltzmannMasslessQuarks


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'

        
def plot_entropy_density_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filename: str,
    parameter_set_annotation: str,
    plotname: str
) -> None:
    print("Building plot: entropy density versus temperature.")
    print(f"Using datafile {filename}.\n")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    # Plot data
    data = FixedChemPotTempData(path_data_folder + filename)
    
    color = 'black'
    linestyle = '-'

    ax.plot(
        data.get_temperature(), 
        data.get_entropy_density(), 
        label=parameter_set_annotation, 
        color=color, 
        linewidth=2, 
        linestyle=linestyle
    )

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    xmin = 0.000
    xmax = 0.300
    ymin = 0.0
    ymax = 0.4
    x_num_ticks = 4
    y_num_ticks = 5
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

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Add text annotations
    auxH = 0.06; auxX = 0.68; auxY = 0.05
    texts = [
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ]
    add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

    fig.tight_layout()

    plt.savefig(path_plots_folder + plotname)

    # Clean up
    plt.clf()
    plt.close()


def plot_entropy_density_vs_temperature_given_list(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filenames: list[str],
    parameter_sets_annotation: list[str],
    colors: list[str],
    linestyles: list[str],
    plotname: str
) -> None:
    print("Building plot: comparing entropy density versus temperature for different files.\n")

    # Check all list sizes
    lengths = [
        len(filenames), 
        len(parameter_sets_annotation),
        len(colors),
        len(linestyles)
    ]
    if (min(lengths)!=max(lengths)):
        print("Lists provided have different sizes!")
        return
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    # Plot data
    for i in range(0, len(filenames)):
        data = FixedChemPotTempData(path_data_folder + filenames[i])

        ax.plot(
            data.get_temperature(), 
            data.get_entropy_density(), 
            label=fr'{parameter_sets_annotation[i]}', 
            color=colors[i], 
            linewidth=2, 
            linestyle=linestyles[i]
        )
    
    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    xmin = 0.1
    xmax = 0.2
    ymin = 0.0
    ymax = 0.105
    x_num_ticks = 5
    y_num_ticks = 6
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
        tick_length=6,
        y_tick_max=0.1
    )

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # Add text annotations
    auxH = 0.06; auxX = 0.68; auxY = 0.05
    texts = [
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ]
    add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

    fig.tight_layout()

    plt.savefig(path_plots_folder + plotname)

    # Clean up
    plt.clf()
    plt.close()

def plot_entropy_density_dPdT_vs_temperature(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filename: str,
    parameter_set_annotation: str,
    plotname: str
) -> None:
    print("Building plot: entropy density, dPdT versus temperature.")
    print(f"Using datafile {filename}.\n")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    # Plot data
    data = FixedChemPotTempData(path_data_folder + filename)
    
    color_entropy = 'black'
    linestyle_entropy = '-'

    color_dPdT = 'red'
    linestyle_dPdT = '--'

    dPdT = np.gradient(data.get_pressure(), data.get_temperature())

    ax.plot(
        data.get_temperature(), 
        data.get_entropy_density(), 
        label=r'$s$', 
        color=color_entropy, 
        linewidth=2, 
        linestyle=linestyle_entropy
    )

    ax.plot(
        data.get_temperature(), 
        dPdT, 
        label=r'$ ({\partial P}/{\partial T})|_{\mu}$', 
        color=color_dPdT, 
        linewidth=2, 
        linestyle=linestyle_dPdT
    )

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$s\, [\mathrm{GeV}^3]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    xmin = 0.000
    xmax = 0.300
    ymin = 0.0
    ymax = 0.4
    x_num_ticks = 4
    y_num_ticks = 5
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

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Add text annotations
    auxH = 0.06; auxX = 0.68; auxY = 0.05
    texts = [
        parameter_set_annotation,
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ]
    add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX, auxY, auxH, texts=texts, fontsize=16)

    fig.tight_layout()

    plt.savefig(path_plots_folder + plotname)

    # Clean up
    plt.clf()
    plt.close()


def plot_s_over_temp3_vs_temp(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    data_specs: list[tuple[str, str, str, int, str]],
    path_output_plot: str,
    stefan_boltzmann_limit: bool = True,
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
        sb = StefanBoltzmannMasslessQuarks(number_of_colors=3,number_of_flavors=3)
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
        ax.legend(loc=legend_loc, fontsize=16, frameon=False)

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$s / T^3$', fontsize=20)
    
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
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
    "set A",
    "entropy_vs_temp_CP0_setA.png"
)

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
    "set B",
    "entropy_vs_temp_CP0_setB.png"
)

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    "set C",
    "entropy_vs_temp_CP0_setC.png"
)

plot_entropy_density_vs_temperature_given_list(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    [
        "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
        "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
        "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    ],
    [
        "set A",
        "set B",
        "set C",
    ],
    [
        'black', 
        'red', 
        'blue',
        ],
    [
        '-', 
        '-', 
        '-',
    ],
    "entropy_vs_temp_CP0_setsABC.png"
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
    "set A",
    "entropy_dPdT_vs_temp_CP0_setA.png"
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
    "set B",
    "entropy_dPdT_vs_temp_CP0_setB.png"
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    "set C",
    "entropy_dPdT_vs_temp_CP0_setC.png"
)

datasets = [
    (
        path_data_folder + "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",  
        r"set A", 
        "black", 
        2, 
        "-"
    ),
    (
        path_data_folder + "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",  
        r"set B", 
        "red", 
        2, 
        "-"
    ),
        (
        path_data_folder + "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat", 
        r"set C", 
        "blue", 
        2, 
        "-"
    ),
]

plot_s_over_temp3_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets,
    path_plots_folder + "s_over_temp3_vs_temp_CP0_setsABC.png",
    stefan_boltzmann_limit=True,
    legend_loc="lower right",
    xlim=(0.120, 0.300),
    ylim=(6.0, 15.0),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.2f",
    y_formatter="%.0f",
    annotation_texts=[
        r"$\mu = 0.0\ \mathrm{GeV}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05
)
