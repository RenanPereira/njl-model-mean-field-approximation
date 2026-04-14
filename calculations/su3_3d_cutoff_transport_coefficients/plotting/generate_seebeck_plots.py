from su3_3d_cutoff_transport_coefficients.plotting.seebeck_plots import (
    plot_seebeck_vs_temp
)


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"


datasets_cpcep = [
    (
        path_transport_data_folder + "SeebeckSigmaeProduct_setA_COMPLETE_COV_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CPCEP.dat", 
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "SeebeckSigmaeProduct_setA_KLEVANSKY_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_KLEVANSKY_CPCEP.dat", 
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "SeebeckSigmaeProduct_setA_ZHUANG_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_ZHUANG_CPCEP.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

plot_seebeck_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "seebeck_vs_temp_setA_CPCEP.png",
    "upper right",
    xlim=(0.040, 0.300),
    ylim=(-9.0, 0.0),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.59,
    y_annotation=0.665,
)
