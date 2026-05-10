from su3_3d_cutoff_quark_relaxation_times.plotting.quark_rel_time_plots import plot_quark_rel_time_vs_temperature
from su3_3d_cutoff_klev_1996_reproduction.plotting.plot_image_overlay import overlay_figure_on_image, plot_quark_rel_time_vs_temperature_over_image


fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
path_data_folder = "su3_3d_cutoff_klev_1996_reproduction/data/"
path_plots_folder = "su3_3d_cutoff_klev_1996_reproduction/plots/"


# Plot quark relaxation time
fig, ax = plot_quark_rel_time_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CP0_B0KlevanskyRecipe.dat",
            "up_quark",
            r"$\tau_{l,\overline{l}}$  | Method II", 
            "black", 
            2, 
            "-"
        ),
        (
            path_data_folder + f"RelaxationTimes_setA_KLEVANSKY_CP0_B0KlevanskyRecipe.dat",
            "strange_quark",
            r"$\tau_{s,\overline{s}}$ | Method II", 
            "red", 
            2, 
            "-"
        ),
    ],
    path_output_plot=path_plots_folder + f"quarks_rel_time_setA_CP0_B0KlevanskyRecipe.png",
    legend_loc="upper right",
    xlim=(0.150, 0.250),
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

# Plot quark relaxation time overlaid with results from paper

background_img_filepath = "su3_3d_cutoff_klev_1996_reproduction/plots/Klev1996NuclPhysAFig15.png"
path_output_image = path_plots_folder + f"quarks_rel_time_setA_CP0_B0KlevanskyRecipe_comparison_simple.png"
path_output_plot = path_plots_folder + f"quarks_rel_time_setA_CP0_B0KlevanskyRecipe_comparison.png"

composed_image = overlay_figure_on_image(
    fig, 
    ax, 
    300, 
    background_img_filepath, 
    (500, 500),
    0.7, 
    path_output_image
)

plot_quark_rel_time_vs_temperature_over_image(
    composed_image, 
    path_output_plot,
    xlim=(0.150, 0.250),
    ylim=(0.0, 10.0),
    x_num_ticks=11,
    y_num_ticks=11,
    tick_fontsize=12,
    x_formatter="%.2f",
    y_formatter="%.0f",
)
