import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import configure_axes, add_annotation_block
from common_utils.shear_viscosity_data import ShearViscosityData
from common_utils.su3_njl_3d_cutoff_data import FixedChemPotTempData
from scipy.interpolate import interp1d
import numpy as np
import math


# Common configurations between plots
#Select font that will be used for the different plots
plt.rcParams['font.family'] = 'sans-serif'


def plot_shear_viscosity_vs_temperature_diff_methods(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filename_complete: str,
    filename_klevansky: str,
    filename_zhuang: str,
    parameter_set_annotation: str,
    plotname: str
) -> None:
    print("Building plot: shear viscosity versus temperature for different methods for the integrated cross section.")
    print(f"Using datafile complete method: {filename_complete}")
    print(f"Using datafile Klevansky method: {filename_klevansky}")
    print(f"Using datafile Zhuang method: {filename_zhuang}\n")
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    # Shear viscosity data
    data_complete = ShearViscosityData(path_data_folder + filename_complete)
    data_klevansky = ShearViscosityData(path_data_folder + filename_klevansky)
    data_zhuang = ShearViscosityData(path_data_folder + filename_zhuang)

    color_complete = 'black'
    linestyle_complete = '-'

    color_klevansky = 'red'
    linestyle_klevansky = '-'

    color_zhuang = 'blue'
    linestyle_zhuang = '-'

    ax.plot(
        data_complete.get_temperature(), 
        data_complete.get_shear_viscosity(), 
        label=r'Method I', 
        color=color_complete, 
        linewidth=2, 
        linestyle=linestyle_complete
    )

    ax.plot(
        data_klevansky.get_temperature(), 
        data_klevansky.get_shear_viscosity(), 
        label=r'Method II', 
        color=color_klevansky, 
        linewidth=2, 
        linestyle=linestyle_klevansky
    )

    ax.plot(
        data_zhuang.get_temperature(), 
        data_zhuang.get_shear_viscosity(), 
        label=r'Method III', 
        color=color_zhuang, 
        linewidth=2, 
        linestyle=linestyle_zhuang
    )

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\eta\, [\mathrm{GeV}^3]$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper left", fontsize=16, frameon=False)

    # Configure axes using the helper function
    xmin = 0.120
    xmax = 0.300
    ymin = 0.0
    ymax = 0.4
    x_num_ticks = 4
    y_num_ticks = 5
    configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize=16, spine_width=1.5, tick_width=1.5, tick_length=6)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Add text annotations
    auxH = 0.06; auxX = 0.03; auxY = 0.3
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


def plot_eta_entropy_ratio_vs_temperature_diff_methods(
    fig_dpi: int,
    fig_x_size: int,
    fig_y_size: int,
    path_data_folder: str,
    path_plots_folder: str,
    filename_complete: str,
    filename_klevansky: str,
    filename_zhuang: str,
    path_file_thermodynamics: str,
    parameter_set_annotation: str,
    plotname: str
) -> None:
    print("Building plot: shear viscosity to entropy density ratio versus temperature for different methods for the integrated cross section.")
    print(f"Using datafile complete method: {filename_complete}")
    print(f"Using datafile Klevansky method: {filename_klevansky}")
    print(f"Using datafile Zhuang method: {filename_zhuang}\n")
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(fig_x_size, fig_y_size), dpi=fig_dpi)
    
    # Shear viscosity data
    data_complete = ShearViscosityData(path_data_folder + filename_complete)
    data_klevansky = ShearViscosityData(path_data_folder + filename_klevansky)
    data_zhuang = ShearViscosityData(path_data_folder + filename_zhuang)

    # Get entropy from thermodynamics data (for this parameter set) and interpolate it
    data_thermodynamics = FixedChemPotTempData(path_file_thermodynamics)
    entropy_dens_interpolation = interp1d(
        data_thermodynamics.get_temperature(), 
        data_thermodynamics.get_entropy_density(), 
        kind='linear'
    )
    entropy_dens = entropy_dens_interpolation(data_complete.get_temperature())

    color_complete = 'black'
    linestyle_complete = '-'

    color_klevansky = 'red'
    linestyle_klevansky = '-'

    color_zhuang = 'blue'
    linestyle_zhuang = '-'

    ax.plot(
        data_complete.get_temperature(), 
        data_complete.get_shear_viscosity()/entropy_dens, 
        label=r'Method I', 
        color=color_complete, 
        linewidth=2, 
        linestyle=linestyle_complete
    )

    ax.plot(
        data_klevansky.get_temperature(), 
        data_klevansky.get_shear_viscosity()/entropy_dens, 
        label=r'Method II', 
        color=color_klevansky, 
        linewidth=2, 
        linestyle=linestyle_klevansky
    )

    ax.plot(
        data_zhuang.get_temperature(), 
        data_zhuang.get_shear_viscosity()/entropy_dens, 
        label=r'Method III', 
        color=color_zhuang, 
        linewidth=2, 
        linestyle=linestyle_zhuang
    )
    
    # Kovtun-Son-Starinets conformal limit (KSS), eta/s|KSS = 1/4*pi
    ax.axhline(
        y=1.0/(4.0*math.pi),
        color='black',
        linestyle='--',
        linewidth=2,
        label=r'KSS'
    )

    # Axes labels
    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=20)
    ax.set_ylabel(r'$\eta/s$', fontsize=20)

    # Grid and legend
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc="upper right", fontsize=16, frameon=False)

    # Configure axes using the helper function
    xmin = 0.120
    xmax = 0.300
    ymin = 0.0
    ymax = 3.0
    x_num_ticks = 4
    y_num_ticks = 7
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
    auxH = 0.06; auxX = 0.05; auxY = 0.85
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


##########################################################################


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_plots_folder = "su3_3d_cutoff_transport_coefficients/plots/"

plot_shear_viscosity_vs_temperature_diff_methods(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "ShearViscosity_setA_COMPLETE_COV.dat",
    "ShearViscosity_setA_KLEVANSKY.dat",
    "ShearViscosity_setA_ZHUANG.dat",
    "set A",
    "shear_viscosity_vs_temp_CP0_setA.png"
)

plot_eta_entropy_ratio_vs_temperature_diff_methods(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    path_data_folder,
    path_plots_folder,
    "ShearViscosity_setA_COMPLETE_COV.dat",
    "ShearViscosity_setA_KLEVANSKY.dat",
    "ShearViscosity_setA_ZHUANG.dat",
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p000000_TMax0p500000_CP0.dat",
    "set A",
    "shear_viscosity_ratio_vs_temp_CP0_setA.png"
)
