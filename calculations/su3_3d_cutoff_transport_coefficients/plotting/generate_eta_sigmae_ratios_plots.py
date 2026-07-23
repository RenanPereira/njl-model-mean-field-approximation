import matplotlib.pyplot as plt
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

####################################################################################################
# set A

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

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
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
plt.close(fig)

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
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
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
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
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
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
plt.close(fig)

datasets_cpcep = [
    (
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CPCEP.dat", 
        "", 
        "black", 
        2, 
        "-"
    ),
]

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_SetA_TMin0p0_TMax0p5_CPU0p318436.dat",
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
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
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
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
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
plt.close(fig)

####################################################################################################
# set B

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setB_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(9.0, 18),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)
plt.close(fig)

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setB_CP0_zoom.png",
    "upper right",
    xlim=(0.168, 0.180),
    ylim=(10.6, 10.9),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setB_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(130, 180),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setB_CP0_zoom.png",
    "upper right",
    xlim=(0.160, 0.180),
    ylim=(134, 140),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p23103.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setB_CPCEP.png",
    "upper right",
    xlim=(0.075, 0.300),
    ylim=(4.0, 28),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setB_CPCEP.png",
    "upper right",
    xlim=(0.040, 0.300),
    ylim=(130, 300),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setB_CPCEP_zoom.png",
    "upper right",
    xlim=(0.097, 0.102),
    ylim=(180, 220),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

####################################################################################################
# set C

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setC_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(9.0, 18),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)
plt.close(fig)

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setC_CP0_zoom.png",
    "upper right",
    xlim=(0.148, 0.160),
    ylim=(10.8, 11.4),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setC_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(130, 180),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CP0.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setC_CP0_zoom.png",
    "upper right",
    xlim=(0.140, 0.170),
    ylim=(134, 140),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_temp_over_sigmae_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p164012.dat",
    path_output_plot_folder + "eta_temp_over_sigmae_s_vs_temp_methods_setC_CPCEP.png",
    "upper right",
    xlim=(0.084, 0.300),
    ylim=(4.0, 36),
    x_num_ticks=5,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.035,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setC_CPCEP.png",
    "upper right",
    xlim=(0.084, 0.300),
    ylim=(130, 300),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)

fig, ax = plot_eta_over_sigmae_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_over_sigmae_temp2_vs_temp_methods_setC_CPCEP_zoom.png",
    "upper right",
    xlim=(0.109, 0.120),
    ylim=(150, 210),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
)
plt.close(fig)
