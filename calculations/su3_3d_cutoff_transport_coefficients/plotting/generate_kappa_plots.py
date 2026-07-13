from su3_3d_cutoff_transport_coefficients.plotting.kappa_plots import (
    plot_kappa_vs_temp, 
    plot_kappa_over_temp2_vs_temp,
)


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"

####################################################################################################
# set A

datasets_cpcep = [
    (
        path_transport_data_folder + "ThermalConductivity_setA_COMPLETE_COV_CPCEP.dat",  
        "", 
        "black", 
        2, 
        "-"
    ),
]

plot_kappa_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "kappa_vs_temp_setA_CPCEP.png",
    "upper left",
    xlim=(0.040, 0.300),
    ylim=(0.0, 15),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

plot_kappa_over_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "kappa_over_temp2_vs_temp_setA_CPCEP.png",
    "upper left",
    xlim=(0.040, 0.300),
    ylim=(0.00, 180.0),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

####################################################################################################
# set B

plot_kappa_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setB_COMPLETE_COV_CPCEP.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_vs_temp_setB_CPCEP.png",
    "upper left",
    xlim=(0.075, 0.300),
    ylim=(0.0, 25),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

plot_kappa_over_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setB_COMPLETE_COV_CPCEP.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_over_temp2_vs_temp_setB_CPCEP.png",
    "upper left",
    xlim=(0.075, 0.300),
    ylim=(0.00, 300.0),
    x_num_ticks=6,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

####################################################################################################
# set C

plot_kappa_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setC_COMPLETE_COV_CPCEP.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_vs_temp_setC_CPCEP.png",
    "upper left",
    xlim=(0.084, 0.300),
    ylim=(0.0, 40),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

plot_kappa_over_temp2_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setC_COMPLETE_COV_CPCEP.dat",  
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_over_temp2_vs_temp_setC_CPCEP.png",
    "upper left",
    xlim=(0.084, 0.300),
    ylim=(0.00, 460.0),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)
