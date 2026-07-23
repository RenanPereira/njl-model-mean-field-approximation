import matplotlib.pyplot as plt
from su3_3d_cutoff_transport_coefficients.plotting.seebeck_plots import (
    plot_seebeck_vs_temp
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
        path_transport_data_folder + "SeebeckSigmaeProduct_setA_COMPLETE_COV_CPCEP.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CPCEP.dat", 
        "", 
        "black", 
        2, 
        "-"
    ),
]

fig, ax = plot_seebeck_vs_temp(
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
plt.close(fig)

####################################################################################################
# set B

fig, ax = plot_seebeck_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "SeebeckSigmaeProduct_setB_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "seebeck_vs_temp_setB_CPCEP.png",
    "upper right",
    xlim=(0.075, 0.300),
    ylim=(-12.0, 0.0),
    x_num_ticks=6,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.59,
    y_annotation=0.665,
)
plt.close(fig)

####################################################################################################
# set C

fig, ax = plot_seebeck_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "SeebeckSigmaeProduct_setC_COMPLETE_COV_CPCEP.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CPCEP.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "seebeck_vs_temp_setC_CPCEP.png",
    "upper right",
    xlim=(0.084, 0.300),
    ylim=(-15.0, 0.0),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.59,
    y_annotation=0.665,
)
plt.close(fig)
