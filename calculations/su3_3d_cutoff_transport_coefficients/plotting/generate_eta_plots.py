from su3_3d_cutoff_transport_coefficients.plotting.eta_plots import plot_eta_vs_temp, plot_eta_over_s_vs_temp


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_transport_data_folder = "su3_3d_cutoff_transport_coefficients/data/"
path_output_plot_folder = "su3_3d_cutoff_transport_coefficients/plots/"

datasets = [
    (
        path_transport_data_folder + "ShearViscosity_setA_COMPLETE_COV_CP0.dat",  
        r"Method I", 
        "black", 
        2, 
        "-"
    ),
    (
        path_transport_data_folder + "ShearViscosity_setA_KLEVANSKY_CP0.dat",  
        r"Method II", 
        "red", 
        2, 
        "-"
    ),
        (
        path_transport_data_folder + "ShearViscosity_setA_ZHUANG_CP0.dat", 
        r"Method III", 
        "blue", 
        2, 
        "-"
    ),
]

path_file_thermo = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_setA_TMin0p000000_TMax0p500000_CP0.dat"

plot_eta_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets,
    path_output_plot_folder + "eta_vs_temp_setA_CP0.png",
    "upper left",
    xlim=(0.120, 0.300),
    ylim=(0.0, 0.4),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu = 0.0\ \mathrm{GeV}$",
    ],
    x_annotation=0.03,
    y_annotation=0.30,
)

plot_eta_over_s_vs_temp(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    datasets,
    path_file_thermo,
    path_output_plot_folder + "eta_over_s_vs_temp_setA_CP0.png",
    include_kss_bound=True,
    legend_loc="upper right",
    xlim=(0.120, 0.300),
    ylim=(0.0, 3.0),
    x_num_ticks=4,
    y_num_ticks=7,
    x_formatter="%.2f", 
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu = 0.0\ \mathrm{GeV}$",
    ],
    x_annotation=0.05,
    y_annotation=0.85,
)
