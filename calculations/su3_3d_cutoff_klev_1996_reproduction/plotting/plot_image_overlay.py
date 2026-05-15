
import matplotlib.pyplot as plt
from PIL import Image
import io
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from common_utils.plot_helper import configure_axes
from matplotlib.figure import Figure
from matplotlib.axes import Axes


# Select font that will be used for the different plots
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False
})


def overlay_figure_on_image(
    fig: Figure, 
    ax: Axes, 
    fig_dpi: int = 150,
    background_img_filepath: str = "",
    size: tuple[float, float] = (500 , 500),
    transparency_factor: float = 0.7,
    path_output_image: str = ""
) -> np.ndarray:
    # Remove ticks and labels from plot
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Convert figure to image (in memory)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=fig_dpi, bbox_inches='tight', pad_inches=0)
    buf.seek(0)
    plot_img = Image.open(buf).convert("RGBA")

    # Load background PNG
    bg_img = Image.open(background_img_filepath).convert("RGBA")

    plot_img = plot_img.resize(size)
    bg_img = bg_img.resize(size)

    # Apply transparency to foreground (your plot)
    plot_img.putalpha(int(255 * transparency_factor))

    # Compose images (overlay)
    composed = Image.alpha_composite(bg_img, plot_img)

    # Save
    composed.save(path_output_image)

    return np.array(composed) 


def plot_quark_rel_time_vs_temperature_over_image(
    composed_image: np.ndarray,
    path_output_plot: str,
    xlim: tuple[float, float] = (0.0 , 1.0),
    ylim: tuple[float, float] = (0.0 , 1.0),
    labels_fontsize: int = 22,
    x_num_ticks: int = 11,
    y_num_ticks: int = 11,
    tick_fontsize: int = 12,
    x_formatter: str = "%.2f",
    y_formatter: str = "%.0f"
) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()

    ax.imshow(
        composed_image,
        extent=[
            xlim[0], 
            xlim[1], 
            ylim[0], 
            ylim[1]
        ],
        origin='upper',
        aspect='auto'
    )

    # Configure axes using the helper function
    xmin = xlim[0]
    xmax = xlim[1]
    ymin = ylim[0]
    ymax = ylim[1]
    configure_axes(
        ax, 
        xmin, 
        xmax, 
        ymin, 
        ymax, 
        x_num_ticks, 
        y_num_ticks, 
        tick_fontsize=tick_fontsize, 
        spine_width=1.5, 
        tick_width=1.5, 
        tick_length=6
    )
    ax.xaxis.set_major_formatter(FormatStrFormatter(x_formatter))
    ax.yaxis.set_major_formatter(FormatStrFormatter(y_formatter))

    ax.set_xlabel(r'$T\, [\mathrm{GeV}]$', fontsize=labels_fontsize)
    ax.set_ylabel(r'$\tau\, [\mathrm{fm}]$', fontsize=labels_fontsize)

    fig.tight_layout()

    fig.savefig(path_output_plot)
    
    return fig, ax
