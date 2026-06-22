from su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.energy_density import (
    plot_energy_vs_temperature, 
    plot_energy_euler_eq_vs_temperature,
    plot_energy_vs_baryon_chem_pot
)

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"

plot_energy_vs_temperature(
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
    path_plots_folder + "energy_vs_temp_setA_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_vs_temperature(
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
    path_plots_folder + "energy_vs_temp_setB_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_vs_temperature(
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
    path_plots_folder + "energy_vs_temp_setC_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_vs_temperature(
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
    path_plots_folder + "energy_vs_temp_setsABC_CP0.png",
    "upper left",
    xlim=(0.1, 0.2),
    ylim=(0.0, 0.015),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.3f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_euler_eq_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon = -P + T \, s + \sum_i \, \mu_i \rho_i$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "energy_euler_eq_vs_temp_setA_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_euler_eq_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon = -P + T \, s + \sum_i \, \mu_i \rho_i$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "energy_euler_eq_vs_temp_setB_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set B",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_euler_eq_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon$', 
            "black", 
            2, 
            "-",
            False
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            r'$\epsilon = -P + T \, s + \sum_i \, \mu_i \rho_i$', 
            "red", 
            2, 
            "--",
            True
        ),
    ],
    path_plots_folder + "energy_euler_eq_vs_temp_setC_CP0.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set C",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.05,
)

plot_energy_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            r'set A , $T\, [\mathrm{GeV}] = 0.075$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetB_T0p12_CPUMin0p0_CPUMax0p5.dat",
            r'set B , $T\, [\mathrm{GeV}] = 0.120$', 
            "red", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetC_T0p13_CPUMin0p0_CPUMax0p5.dat",
            r'set C , $T\, [\mathrm{GeV}] = 0.130$', 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "energy_vs_muB_setsABC_diffT.png",
    "upper left",
    xlim=(0.0, 1.500),
    ylim=(0.0, 0.036),
    x_num_ticks=7,
    y_num_ticks=7,
    x_formatter="%.2f", 
    y_formatter="%.3f",
)
