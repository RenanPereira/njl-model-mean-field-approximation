from su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.entropy_density import (
    plot_entropy_density_vs_temperature, 
    plot_entropy_density_dPdT_vs_temperature,
    plot_s_over_temp3_vs_temp
)

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
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            "set A", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "entropy_vs_temp_CP0_setA.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            "set B", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "entropy_vs_temp_CP0_setB.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            "set C", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "entropy_vs_temp_CP0_setC.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            "set A", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            "set B", 
            "red", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            "set C", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "entropy_vs_temp_CP0_setsABC.png",
    "upper left",
    xlim=(0.1, 0.2),
    ylim=(0.0, 0.105),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.2f", 
    y_formatter="%.2f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$s$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$ ({\partial P}/{\partial T})|_{\mu}$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "entropy_dPdT_vs_temp_CP0_setA.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$s$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$ ({\partial P}/{\partial T})|_{\mu}$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "entropy_dPdT_vs_temp_CP0_setB.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_entropy_density_dPdT_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$s$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$ ({\partial P}/{\partial T})|_{\mu}$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "entropy_dPdT_vs_temp_CP0_setC.png",
    "upper left",
    xlim=(0.0, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_s_over_temp3_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
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
    ],
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
