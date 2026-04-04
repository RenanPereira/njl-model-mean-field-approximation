from su3_3d_cutoff_transport_coefficients.plotting.sigmae_plots import (
    plot_sigmae_vs_temp, 
    plot_sigmae_over_temp_vs_temp
)


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"

datasets_cp0 = [
    (
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CP0.dat",  
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY_CP0.dat",  
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG_CP0.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_sigmae_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    path_output_plot_folder + "sigmae_vs_temp_setA_CP0.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.03),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

plot_sigmae_over_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cp0,
    path_output_plot_folder + "sigmae_over_temp_vs_temp_setA_CP0.png",
    "upper right",
    xlim=(0.120, 0.300),
    ylim=(0.00, 0.15),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.03,
    y_annotation=0.05,
)

datasets_cpcep = [
    (
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CPCEP.dat",  
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY_CPCEP.dat",  
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG_CPCEP.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_sigmae_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "sigmae_vs_temp_setA_CPCEP.png",
    "upper left",
    xlim=(0.040, 0.300),
    ylim=(0.0, 0.1),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)

plot_sigmae_over_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "sigmae_over_temp_vs_temp_setA_CPCEP.png",
    "upper right",
    xlim=(0.040, 0.300),
    ylim=(0.00, 0.9),
    x_num_ticks=5,
    y_num_ticks=6,
    x_formatter="%.3f", 
    y_formatter="%.2f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.03,
    y_annotation=0.05,
)
