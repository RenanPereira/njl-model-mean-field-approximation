import matplotlib.pyplot as plt
from su3_3d_cutoff_transport_coefficients.plotting.eta_plots import plot_eta_vs_temp, plot_eta_over_s_vs_temp


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
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CP0p0.dat",  
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ShearViscosity_setA_KLEVANSKY_CP0p0.dat",  
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ShearViscosity_setA_ZHUANG_CP0p0.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    path_output_plot_folder + "eta_vs_temp_setA_CP0p0.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setA_CP0p0.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.120, 0.300),
    ylim=(0.0, 3.0),
    x_num_ticks=4,
    y_num_ticks=7,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)

datasets_cpcep = [
    (
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CP0p318436.dat",  
        "", 
        "black", 
        2, 
        "-"
    ),
]

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "eta_vs_temp_setA_CP0p318436.png",
    "upper left",
    xlim=(0.040, 0.300),
    ylim=(0.0, 1.5),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p0_TMax0p5_CPU0p318436.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setA_CP0p318436.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.040, 0.300),
    ylim=(0.0, 6.0),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)

####################################################################################################
# set B

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0p0.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_vs_temp_setB_CP0p0.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0p0.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setB_CP0p0.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.120, 0.300),
    ylim=(0.0, 3.0),
    x_num_ticks=4,
    y_num_ticks=7,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0p231030.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_vs_temp_setB_CP0p231030.png",
    "upper left",
    xlim=(0.075, 0.300),
    ylim=(0.0, 1.5),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setB_COMPLETE_COV_CP0p231030.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setB_TMin0p0_TMax0p5_CPU0p23103.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setB_CP0p231030.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.075, 0.300),
    ylim=(0.0, 6.0),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)

####################################################################################################
# set C

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0p0.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_vs_temp_setC_CP0p0.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0p0.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p0.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setC_CP0p0.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.120, 0.300),
    ylim=(0.0, 3.0),
    x_num_ticks=4,
    y_num_ticks=7,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)

fig, ax = plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0p164012.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "eta_vs_temp_setC_CP0p164012.png",
    "upper left",
    xlim=(0.084, 0.300),
    ylim=(0.0, 1.5),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.60,
    y_annotation=0.88,
)
plt.close(fig)

fig, ax = plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ShearViscosity_setC_COMPLETE_COV_CP0p164012.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setC_TMin0p0_TMax0p5_CPU0p164012.dat",
    path_output_plot_folder + "eta_over_s_vs_temp_setC_CP0p164012.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.084, 0.300),
    ylim=(0.0, 6.0),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
plt.close(fig)
