import matplotlib.pyplot as plt
from su3_3d_cutoff_transport_coefficients.plotting.kappa_sigmae_ratio_plots import (
    plot_kappa_over_sigmae_temp_vs_temp
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
        path_transport_data_folder + "ThermalConductivity_setA_COMPLETE_COV_CP0p318436.dat", 
        path_transport_data_folder + "ElectricalConductivity_setA_COMPLETE_COV_CP0p318436.dat", 
        "", 
        "black", 
        2, 
        "-"
    ),
]

fig, ax = plot_kappa_over_sigmae_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets_cpcep,
    path_output_plot_folder + "kappa_over_sigmae_temp_vs_temp_setA_CP0p318436.png",
    "upper left",
    xlim=(0.040, 0.300),
    ylim=(100, 600),
    x_num_ticks=5,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)
plt.close(fig)

####################################################################################################
# set B

fig, ax = plot_kappa_over_sigmae_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setB_COMPLETE_COV_CP0p231030.dat", 
            path_transport_data_folder + "ElectricalConductivity_setB_COMPLETE_COV_CP0p231030.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_over_sigmae_temp_vs_temp_setB_CP0p231030.png",
    "upper left",
    xlim=(0.075, 0.300),
    ylim=(100, 1100),
    x_num_ticks=6,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set B",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)
plt.close(fig)

####################################################################################################
# set C

fig, ax = plot_kappa_over_sigmae_temp_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    [
        (
            path_transport_data_folder + "ThermalConductivity_setC_COMPLETE_COV_CP0p164012.dat", 
            path_transport_data_folder + "ElectricalConductivity_setC_COMPLETE_COV_CP0p164012.dat", 
            "", 
            "black", 
            2, 
            "-"
        ),
    ],
    path_output_plot_folder + "kappa_over_sigmae_temp_vs_temp_setC_CP0p164012.png",
    "upper left",
    xlim=(0.084, 0.300),
    ylim=(100, 2400),
    x_num_ticks=5,
    y_num_ticks=5,
    x_formatter="%.3f", 
    y_formatter="%.0f",
    annotation_texts=[
        "set C",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.05,
    y_annotation=0.665,
)
plt.close(fig)
