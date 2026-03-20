from su3_3d_cutoff_int_cross_sections.plotting_utils import plot_integrated_cross_section_vs_temperature
from common_utils.cross_section_data import process_to_ylabel_latex

# Common configurations between plots
fig_dpi = 150
fig_x_size = 6
fig_y_size = 6

# Location of the data and plots folder with respect to calculations folder
data_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/data/"
plots_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/plots/"

process = "UUUU"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "UUBarUUBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 12.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "UUBarDDBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "UUBarSSBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 3.0),
    x_num_ticks=4,
    y_num_ticks=4,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.08,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "USUS"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "USBarUSBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "UDUD"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "UDBarUDBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="lower left",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "SSSS"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)

process = "SSBarUUBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 8.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.05,
    annotation_vert_space=0.06
)

process = "SSBarSSBar"
plot_integrated_cross_section_vs_temperature(
    fig_dpi,
    fig_x_size,
    fig_y_size,
    data_specs=[
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_COMPLETE_COV_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method I", 
            "black", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_KLEVANSKY_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method II", 
            "red", 
            2, 
            "-"
        ),
        (
            [
                data_folder + f"IntegratedCrossSection_setA_{process}_ZHUANG_TMin0p120000_TMax0p300000_CPU0p000000.dat",
            ],  
            r"Method III", 
            "blue", 
            2, 
            "-"
        ),
    ],
    path_output_plot=plots_folder + f"integrated_cross_section_setA_{process}_CP0.png",
    legend_loc="upper right",
    label_int_cross_section=process_to_ylabel_latex(process),
    xlim=(0.120, 0.300),
    ylim=(0.0, 12.0),
    x_num_ticks=4,
    y_num_ticks=5,
    x_formatter="%.2f",
    y_formatter="%.1f",
    annotation_texts=[
        "set A",
        r"$\mu [\mathrm{GeV}] = 0.0$",
    ],
    x_annotation=0.05,
    y_annotation=0.88,
    annotation_vert_space=0.06
)
