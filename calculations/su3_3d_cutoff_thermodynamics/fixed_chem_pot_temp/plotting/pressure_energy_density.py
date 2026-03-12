import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_pressure_vs_energy(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filename: str,
    parameter_set_annotation: str,
    plotname: str
) -> None:
    print("Building plot: pressure versus energy density.")
    print(f"Using datafile {filename}.\n")

    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)

    # Plot data
    data = FixedChemPotTempData(path_data_folder + filename)
    
    color = 'black'
    linestyle = '-'

    ax.plot(
        data.get_energy_density(), 
        data.get_pressure(), 
        label=parameter_set_annotation, 
        color=color, 
        linewidth=2, 
        linestyle=linestyle
    )

    # Axes labels
    ax.set_xlabel(r'$\epsilon \, [\mathrm{GeV}^4]$', fontsize=20)
    ax.set_ylabel(r'$P \, [\mathrm{GeV}^4]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    xmin = 0.000
    xmax = 0.6
    ymin = 0.0
    ymax = 0.2
    x_num_ticks = 5
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


def plot_pressure_vs_energy_given_list(
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
    print("Building plot: comparing pressure versus energy density for different files.\n")

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
            data.get_energy_density(), 
            data.get_pressure(), 
            label=fr'{parameter_sets_annotation[i]}', 
            color=colors[i], 
            linewidth=2, 
            linestyle=linestyles[i]
        )
    
    # Axes labels
    ax.set_xlabel(r'$\epsilon \, [\mathrm{GeV}^4]$', fontsize=20)
    ax.set_ylabel(r'$P \, [\mathrm{GeV}^4]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False, title_fontsize=14)

    # Configure axes using the helper function
    xmin = 0.000
    xmax = 0.6
    ymin = 0.0
    ymax = 0.2
    x_num_ticks = 5
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


##########################################################################


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"

plot_pressure_vs_energy(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p000000_TMax0p500000_CP0.dat",
    "set A",
    "pressure_vs_energy_CP0_setA.png"
)

plot_pressure_vs_energy(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p000000_TMax0p500000_CP0.dat",
    "set B",
    "pressure_vs_energy_CP0_setB.png"
)

plot_pressure_vs_energy(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p000000_TMax0p500000_CP0.dat",
    "set C",
    "pressure_vs_energy_CP0_setC.png"
)

plot_pressure_vs_energy_given_list(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    [
        "SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p000000_TMax0p500000_CP0.dat",
        "SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p000000_TMax0p500000_CP0.dat",
        "SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p000000_TMax0p500000_CP0.dat",
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
    "pressure_vs_energy_CP0_setsABC.png"
)
