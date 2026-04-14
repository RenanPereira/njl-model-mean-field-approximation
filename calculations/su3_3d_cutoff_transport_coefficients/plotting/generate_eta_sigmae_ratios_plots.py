from su3_3d_cutoff_transport_coefficients.plotting.eta_sigmae_ratios_plots import (
    plot_eta_temp_over_sigmae_s_vs_temp, 
    plot_eta_over_sigmae_temp2_vs_temp
)


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"

datasets_cp0 = [
    (
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CP0.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CP0.dat", 
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ShearViscosity_setA_KLEVANSKY_CP0.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY_CP0.dat", 
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ShearViscosity_setA_ZHUANG_CP0.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG_CP0.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setA_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(9.0, 18),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)

plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setA_CP0_zoom.png",
    "upper right",
    xlim=(0.200, 0.225),
    ylim=(10.3, 10.8),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)

plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setA_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(130, 180),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)

plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setA_CP0_zoom.png",
    "upper right",
    xlim=(0.200, 0.225),
    ylim=(135, 140),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)

datasets_cpcep = [
    (
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CPCEP.dat", 
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ShearViscosity_setA_KLEVANSKY_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY_CPCEP.dat", 
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ShearViscosity_setA_ZHUANG_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG_CPCEP.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_SetA_TMin0p0_TMax0p5_CPU0p318434.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setA_CPCEP.png",
    "upper right",
    xlim=(0.040, 0.300),
    ylim=(4.0, 12),
    x_num_ticks=5,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)

plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setA_CPCEP.png",
    "upper right",
    xlim=(0.040, 0.300),
    ylim=(150, 300),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)

plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setA_CPCEP_zoom.png",
    "upper right",
    xlim=(0.065, 0.070),
    ylim=(288, 300),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
