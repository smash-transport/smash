#!/usr/bin/env python3

# ===================================================
#
#    Copyright (c) 2024
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
# ===================================================

import matplotlib.pyplot as plt


def make_multi_bar_histogram_plot(data, colors=None, total_width=0.8,
                                  single_width=1, legend=True,
                                  error_thickness=0.8, xticks=None,
                                  ylabel="", ymargin=0.0, annotate=False):
    """Draws a bar plot with multiple bars per data point.

    Parameters
    ----------
    data: dictionary
        A dictionary containing the data we want to plot. Keys are the names
        of the data, the items is a list of (value, error) pairs.

        Example:
        data = {
            "dataset_1": [(1,0.3), (2,0.8), (3,1)],
            "dataset_2": [(4,0.3), (5,0.8), (6,1)],
            "dataset_3": [(7,0.3), (8,0.8), (9,1)],
        }

    colors : array-like, optional
        A list of colors which are used for the bars. If None, the colors
        will be the standard matplotlib color cycle. (default: None)

    total_width : float, optional, default: 0.8
        The width of a bar group. 0.8 means that 80% of the x-axis is covered
        by bars and 20% will be spaces between the bars.

    single_width: float, optional, default: 1
        The relative width of a single bar within a group. 1 means the bars
        will touch each other within a group, values less than 1 will make
        these bars thinner.

    legend: bool, optional, default: True
        If this is set to true, a legend will be added to the axis.

    error_thickness: float, optional, default: 0.8
        Width of the error line and cap in points.

    xticks : array-like, optional
        A list of labels which are used as x-tick for the bars. If None, these
        will be just integers starting from 0. (default: None)

    ylabel: string, optional, default: ""
        Label of the y-axis.

    ymargin: float, optional, default: 0.0
        The padding added to the top limit of the y-axis relative to its
        natural range.

    annotate: bool or array-like, optional, default: False
        If this is set to True, an annotation on each bar (except the first in
        each group) will be placed. This is in percentage the height difference
        w.r.t. the previous bar in the same group. If the value is array-like
        only the bars in the given list will be annotated, e.g. [1,2] will
        annotate the second and third bar in each group. The value 0 is accepted
        but it does not produce an annotation as it does not make sense.
    """

    # Prepare figure with axis
    fig, ax = plt.subplots()

    # Check if colors where provided, otherwise use the default color cycle
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Number of bars per group
    n_bars = len(data)

    # The width of a single bar
    bar_width = total_width / n_bars

    # List containing handles for the drawn bars, used for the legend
    bars = []

    # Iterate over all data
    for i, (name, values) in enumerate(data.items()):
        # The offset in x direction of that bar
        x_offset = (i - n_bars / 2) * bar_width + bar_width / 2

        # Draw a bar for every value of that type
        for x, y_dy in enumerate(values):
            y = y_dy[0]
            dy = y_dy[1]
            bar = ax.bar(x + x_offset, y, width=bar_width * single_width,
                         yerr=dy, capsize=3, color=colors[i % len(colors)],
                         error_kw={"elinewidth": error_thickness,
                                   "capthick": error_thickness})
            # Understand what to annotate and do so if requested
            if annotate == True:
                bars_to_annotate = [i for i, _ in enumerate(data.keys())]
            elif annotate == False:
                bars_to_annotate = []
            else:
                bars_to_annotate = annotate
            if i in bars_to_annotate and i > 0:
                previous_bar_height = list(data.values())[i-1][x][0]
                if y != 0:
                    delta = (y-previous_bar_height)/y*100
                    plt.annotate(f'{delta:.1f}%', xy=(x + x_offset, y+dy),
                                 xytext=(0, 5), textcoords='offset points',
                                 ha='center', va='bottom', rotation=90,
                                 color='red' if delta > 0 else 'green',
                                 fontsize='small', fontweight='bold')

        # Add a handle to the last drawn bar, which we'll need for the legend
        bars.append(bar[0])

    # Draw legend if we need
    if legend:
        ax.legend(bars, data.keys(),
                  ncols=2, bbox_to_anchor=(0.5, 1), loc='lower center')

    # Customize x-ticks if desired
    if xticks is not None:
        plt.xticks(range(len(next(iter(data.values())))), xticks,
                   rotation=30, ha="right", rotation_mode="anchor")

    # Add label to y-axis
    if ylabel != "":
        plt.ylabel(ylabel)

    # Adjust the y-axis margin if desired
    if ymargin != 0.0:
        ax.margins(y=ymargin)

    # Show the plot
    plt.tight_layout()  # prevent x-ticks and/or axis labels to be cut
    plt.show()


if __name__ == "__main__":
    # Usage example:
    data = {
        "a": [(1, 1), (2, 1), (3, 1), (2, 1), (1, 1)],
        "b": [(2, 1), (3, 1), (4, 1), (3, 1), (1, 1)],
        "c": [(3, 1), (2, 1), (1, 1), (4, 1), (2, 1)]
    }

    make_multi_bar_histogram_plot(data, total_width=.8, single_width=.95,
                                  xticks=["one", "two",
                                          "three", "four", "five"],
                                  ylabel="Magic number", ymargin=0.2,
                                  annotate=True)
