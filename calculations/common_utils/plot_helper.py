import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms


def pos(min_val, max_val, factor):
    return min_val + factor * (max_val - min_val)


def configure_axes(ax, xmin, xmax, ymin, ymax, x_num_ticks, y_num_ticks, tick_fontsize, spine_width, tick_width, tick_length):
    """
    Configures the axes of a matplotlib plot.

    Parameters:
    - ax: The axes object to configure.
    - xmin, xmax: The x-axis limits.
    - ymin, ymax: The y-axis limits.
    - x_num_ticks, y_num_ticks: Number of ticks on x and y axes.
    - tick_fontsize: Font size for the tick labels.
    - spine_width: Line width for the axes spines.
    - tick_width: Width of the ticks.
    - tick_length: Length of the ticks.
    """
    # Set x and y limits
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    # Setup x-axis ticks
    x_ticks = np.linspace(xmin, xmax, x_num_ticks)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{tick:.1f}' for tick in x_ticks], fontsize=tick_fontsize)

    # Setup y-axis ticks
    y_ticks = np.linspace(ymin, ymax, y_num_ticks)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{tick:.1f}' for tick in y_ticks], fontsize=tick_fontsize)

    # Make the axis and ticks lines thicker
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(spine_width)
    ax.tick_params(axis='both', width=tick_width, length=tick_length)  # Customize tick size and thickness


def add_annotation(ax, xmin, xmax, ymin, ymax, auxX, auxY, text, fontsize=14):
    """
    Adds a single annotation to the given axis and returns the absolute rectangle coordinates.

    Parameters:
    - ax: The matplotlib axis object.
    - xmin, xmax: x-axis limits.
    - ymin, ymax: y-axis limits.
    - auxX, auxY: Relative positions for the text placement.
    - text: The annotation string.
    - fontsize: Font size for the annotation.

    Returns:
    - absolute_coords: A tuple of the form (x_abs_min, y_abs_min, x_abs_max, y_abs_max)
                       where the coordinates are in absolute data space.
    """
    # Create the text annotation
    text_obj = ax.text(pos(xmin, xmax, auxX), pos(ymin, ymax, auxY), text, fontsize=fontsize)

    # Draw the figure to compute bounding box
    plt.draw()

    # Get the bounding box in display coordinates
    bbox = text_obj.get_window_extent(renderer=ax.figure.canvas.get_renderer())

    # Transform bounding box to data coordinates
    inv = ax.transData.inverted()
    bbox_data = bbox.transformed(inv)

    # Extract absolute rectangle coordinates
    x_abs_min, y_abs_min = bbox_data.xmin, bbox_data.ymin
    x_abs_max, y_abs_max = bbox_data.xmax, bbox_data.ymax
    absolute_coords = (x_abs_min, y_abs_min, x_abs_max, y_abs_max)

    return absolute_coords


def add_annotation_block(ax, xmin, xmax, ymin, ymax, auxX1, auxY1, auxH, texts, fontsize):
    """
    Adds a block of annotations to the given axis by calling `add_annotation` for each text.

    Parameters:
    - ax: The matplotlib axis object.
    - xmin, xmax: x-axis limits.
    - ymin, ymax: y-axis limits.
    - auxX1, auxY1: Relative positions for the starting text placement.
    - auxH: Vertical spacing between annotations.
    - texts: List of annotation strings.
    - fontsize: Font size for the annotations.
    """
    for i, text in enumerate(texts):
        auxY = auxY1 + auxH * (len(texts) - 1 - i)  # Compute relative y-position for each annotation
        add_annotation(ax, xmin, xmax, ymin, ymax, auxX1, auxY, text, fontsize)


def annotate_with_2_lines(ax, xmin, xmax, ymin, ymax, label_line_dist, text, x_start_factor, y_factor,
                          length_line_fraction_of_box, fontsize=16, color1="black",
                          color2="black", width1=2, width2=2, style1="-", style2="--"):
    """
    Adds a custom annotation to a plot and two lines extending from the annotation box.

    Parameters:
        ax (matplotlib.axes.Axes): The axes to draw on.
        xmin, xmax, ymin, ymax (float): Coordinates of the box area.
        text (str): The annotation text.
        x_start_factor (float): Factor to calculate x start for the annotation box.
        y_factor (float): Factor to calculate the y position of the annotation box.
        length_line_fraction_of_box (float): Fraction of box width to calculate line length.
        fontsize (int): Font size for the text.
        color1, color2 (str): Colors for the two lines.
        width1, width2 (float): Line widths for the two lines.
        style1, style2 (str): Line styles for the two lines.
        label_line_dist: Distance between the label and the start of the lines.

    Returns:
        None
    """
    # Adding annotation and determining rectangle coordinates
    rectangle_coords = add_annotation(ax, xmin, xmax, ymin, ymax, x_start_factor, y_factor, text, fontsize)

    # Box properties
    box_height=(rectangle_coords[3]-rectangle_coords[1])
    box_start = -label_line_dist + rectangle_coords[0]
    box_end = label_line_dist + rectangle_coords[2]

    # Adding the first line
    y_pos = rectangle_coords[1] + box_height*(2./3.)
    ax.plot([box_end, pos(box_end, box_start, -length_line_fraction_of_box)], 
            [y_pos, y_pos], 
            color=color1, linewidth=width1, linestyle=style1)

    # Adding the second line
    y_pos = rectangle_coords[1] + box_height*(1./3.)
    ax.plot([box_end, pos(box_end, box_start, -length_line_fraction_of_box)], 
            [y_pos, y_pos], 
            color=color2, linewidth=width2, linestyle=style2)
