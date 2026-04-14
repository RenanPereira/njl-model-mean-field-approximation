from su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.quark_density_plots import (
    plot_quark_dens_vs_chem_pot, 
    plot_baryon_dens_vs_baryon_chem_pot, 
    plot_baryon_dens_dPdmuB_vs_baryon_chem_pot
)

fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/"
path_plots_folder = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/plots/"

plot_quark_dens_vs_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            "up_quark",
            r'$\rho_l$', 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            "strange_quark",
            r'$\rho_s$', 
            "red", 
            2, 
            "-"
        ),
    ],
    "up_quark",
    path_plots_folder + "rho_quarks_vs_mu_setA_T0p075.png",
    "upper left",
    xlim=(0.0, 0.500),
    ylim=(0.0, 2.0),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.1f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$T\, [\mathrm{GeV}] = 0.075$",
    ],
    x_annotation=0.05,
    y_annotation=0.08,
)

plot_baryon_dens_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_plots_folder + "rhoB_vs_muB_setA_T0p075.png",
    None,
    xlim=(0.0, 1.500),
    ylim=(0.0, 2.0),
    x_num_ticks=7,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$T\, [\mathrm{GeV}] = 0.075$",
    ],
    x_annotation=0.05,
    y_annotation=0.08,
)

plot_baryon_dens_dPdmuB_vs_baryon_chem_pot(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_data_folder + f"SU3NJL3DCutoffFixedChemPotTemp_SetA_T0p075_CPUMin0p0_CPUMax0p5.dat",
            r'$\rho_B$', 
            "black", 
            2, 
            "-",
            r'$ ({\partial P}/{\partial \mu_B})|_{T}$', 
            "red", 
            2, 
            "--"
        ),
    ],
    path_plots_folder + "rhoB_dPdmuB_vs_muB_setA_T0p075.png",
    "upper left",
    xlim=(0.0, 1.500),
    ylim=(0.0, 2.0),
    x_num_ticks=7,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$T\, [\mathrm{GeV}] = 0.075$",
    ],
    x_annotation=0.05,
    y_annotation=0.08,
)
