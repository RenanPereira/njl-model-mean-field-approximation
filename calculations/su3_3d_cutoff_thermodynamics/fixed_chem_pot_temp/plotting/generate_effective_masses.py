from su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.effective_masses import (
    plot_quark_masses_vs_temperature, 
    plot_normalized_quark_masses_vs_temperature,
    plot_quark_masses_vs_baryon_chem_pot
)

# Plots properties
fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"

plot_quark_masses_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_temp_CP0_setA.png",
    "upper right",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=7,
    x_formatter="%.1f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.1,
)

plot_quark_masses_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_temp_CP0_setB.png",
    "upper right",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=7,
    x_formatter="%.1f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.1,
)

plot_quark_masses_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_temp_CP0_setC.png",
    "upper right",
    xlim=(0.0, 0.500),
    ylim=(0.0, 0.6),
    x_num_ticks=6,
    y_num_ticks=7,
    x_formatter="%.1f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.1,
)

plot_normalized_quark_masses_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'set A', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'set B', 
            "red", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
            "up_quark",
            r'set C', 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_temp_CP0_setsABC.png",
    "upper right",
    xlim=(0.0, 0.500),
    ylim=(0.0, 1.05),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.1f", 
    y_formatter="%.1f",
    annotation_texts=[
        r'$\mu = 0.0\ \mathrm{GeV}$',
    ],
    x_annotation=0.68,
    y_annotation=0.1,
)

plot_quark_masses_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_muB_setA_T0p075.png",
    "upper right",
    xlim=(0.0, 1.500),
    ylim=(0.0, 0.6),
    x_num_ticks=7,
    y_num_ticks=7,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$T\, [\mathrm{GeV}] = 0.075$",
    ],
    x_annotation=0.05,
    y_annotation=0.1,
)

plot_quark_masses_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetB_T0p12_CPUMin0p0_CPUMax0p5.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetB_T0p12_CPUMin0p0_CPUMax0p5.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_muB_setB_T0p12.png",
    "upper right",
    xlim=(0.0, 1.500),
    ylim=(0.0, 0.6),
    x_num_ticks=7,
    y_num_ticks=7,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$T\, [\mathrm{GeV}] = 0.120$",
    ],
    x_annotation=0.6,
    y_annotation=0.1,
)

plot_quark_masses_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetC_T0p13_CPUMin0p0_CPUMax0p5.dat",
            "up_quark",
            r'$M_\ell$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetC_T0p13_CPUMin0p0_CPUMax0p5.dat",
            "strange_quark",
            r'$M_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "quark_eff_masses_vs_muB_setC_T0p13.png",
    "upper right",
    xlim=(0.0, 1.500),
    ylim=(0.0, 0.6),
    x_num_ticks=7,
    y_num_ticks=7,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$T\, [\mathrm{GeV}] = 0.130$",
    ],
    x_annotation=0.6,
    y_annotation=0.1,
)
