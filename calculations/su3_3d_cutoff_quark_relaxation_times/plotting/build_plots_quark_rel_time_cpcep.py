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
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "up_quark",
            r"$\tau_{l}$",
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "strange_quark",
            r"$\tau_{s}$",
            "red", 
            2, 
            ":"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "up_antiquark",
            r"$\tau_{\overline{l}}$",
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "strange_antiquark",
            r"$\tau_{\overline{s}}$",
            "red", 
            2, 
            "-."
        ),
    ],
    path_output_plot=path_plots_folder + f"quarks_rel_time_setA_methodI_CPCEP.png",
    legend_loc="upper right",
    xlim=(0.040, 0.300),
    ylim=(0.0, 30.0),
    x_num_ticks=5,
    y_num_ticks=7,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        "Method I",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.10,
    y_annotation=0.8,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "up_quark",
            r"$\tau_{l}$",
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "strange_quark",
            r"$\tau_{s}$",
            "red", 
            2, 
            ":"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "up_antiquark",
            r"$\tau_{\overline{l}}$",
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "strange_antiquark",
            r"$\tau_{\overline{s}}$",
            "red", 
            2, 
            "-."
        ),
    ],
    path_output_plot=path_plots_folder + f"quarks_rel_time_setA_methodII_CPCEP.png",
    legend_loc="upper right",
    xlim=(0.040, 0.300),
    ylim=(0.0, 30.0),
    x_num_ticks=5,
    y_num_ticks=7,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        "Method II",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.10,
    y_annotation=0.8,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "up_quark",
            r"$\tau_{l}$",
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "strange_quark",
            r"$\tau_{s}$",
            "red", 
            2, 
            ":"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "up_antiquark",
            r"$\tau_{\overline{l}}$",
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "strange_antiquark",
            r"$\tau_{\overline{s}}$",
            "red", 
            2, 
            "-."
        ),
    ],
    path_output_plot=path_plots_folder + f"quarks_rel_time_setA_methodIII_CPCEP.png",
    legend_loc="upper right",
    xlim=(0.040, 0.300),
    ylim=(0.0, 30.0),
    x_num_ticks=5,
    y_num_ticks=7,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        "Method III",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.10,
    y_annotation=0.8,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "up_quark",
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "up_quark",
            r"Method II", 
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "up_quark",
            r"Method III", 
            "black", 
            2, 
            ":"
        ),
    ],
    path_output_plot=path_plots_folder + f"l_quark_rel_time_setA_CPCEP.png",
    legend_loc="upper right",
    label_rel_time=r"$\tau_{l} \, [\mathrm{fm}]$",
    xlim=(0.040, 0.300),
    ylim=(0.0, 30.0),
    x_num_ticks=5,
    y_num_ticks=7,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.14,
    y_annotation=0.8+0.06,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "up_antiquark",
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "up_antiquark",
            r"Method II", 
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "up_antiquark",
            r"Method III", 
            "black", 
            2, 
            ":"
        ),
    ],
    path_output_plot=path_plots_folder + f"l_antiquark_rel_time_setA_CPCEP.png",
    legend_loc="upper right",
    label_rel_time=r"$\tau_{\overline{l}} \, [\mathrm{fm}]$",
    xlim=(0.040, 0.300),
    ylim=(0.0, 7.5),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.14,
    y_annotation=0.8+0.06,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "strange_quark",
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "strange_quark",
            r"Method II", 
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "strange_quark",
            r"Method III", 
            "black", 
            2, 
            ":"
        ),
    ],
    path_output_plot=path_plots_folder + f"s_quark_rel_time_setA_CPCEP.png",
    legend_loc="upper right",
    label_rel_time=r"$\tau_{s} \, [\mathrm{fm}]$",
    xlim=(0.040, 0.300),
    ylim=(0.0, 30.0),
    x_num_ticks=5,
    y_num_ticks=7,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.14,
    y_annotation=0.8+0.06,
    annotation_vert_space=0.06
)

plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_COMPLETE_COV_CPCEP.dat",
            "strange_antiquark",
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CPCEP.dat",
            "strange_antiquark",
            r"Method II", 
            "black", 
            2, 
            "--"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_ZHUANG_CPCEP.dat",
            "strange_antiquark",
            r"Method III", 
            "black", 
            2, 
            ":"
        ),
    ],
    path_output_plot=path_plots_folder + f"s_antiquark_rel_time_setA_CPCEP.png",
    legend_loc="upper right",
    label_rel_time=r"$\tau_{\overline{s}} \, [\mathrm{fm}]$",
    xlim=(0.040, 0.300),
    ylim=(0.0, 7.5),
    x_num_ticks=5,
    y_num_ticks=4,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = \mu_{\mathrm{CEP}}$",
    ],
    x_annotation=0.14,
    y_annotation=0.8+0.06,
    annotation_vert_space=0.06
)
