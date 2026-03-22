from su3_3d_cutoff_quark_relaxation_times.plotting.quark_rel_time_plots import plot_quark_rel_time_vs_temperature


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_plots_folder = "su3_3d_cutoff_quark_relaxation_times/plots/"


plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CP0.dat",
            "up_quark",
            r"$\tau_{l,\overline{l}}$  | Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CP0.dat",
            "strange_quark",
            r"$\tau_{s,\overline{s}}$ | Method I", 
            "red", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CP0.dat",
            "up_quark",
            r"$\tau_{l,\overline{l}}$  | Method II", 
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CP0.dat",
            "strange_quark",
            r"$\tau_{s,\overline{s}}$ | Method II", 
            "red", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CP0.dat",
            "up_quark",
            r"$\tau_{l,\overline{l}}$  | Method III", 
            "black", 
            2, 
            ":"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CP0.dat",
            "strange_quark",
            r"$\tau_{s,\overline{s}}$ | Method III", 
            "red", 
            2, 
            ":"
        ),
    ],
    path_output_plot=path_plots_folder + f"quark_relaxation_time_setA_CP0.png",
    legend_loc="upper right",
    xlim=(0.120, 0.300),
    ylim=(0.0, 10.0),
    x_num_ticks=4,
    y_num_ticks=6,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.03,
    y_annotation=0.03,
    annotation_vert_space=0.06
)
