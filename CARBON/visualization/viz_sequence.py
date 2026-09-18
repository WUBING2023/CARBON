"""
Sequence logo and score track visualization (tfmodisco-style).
Adapted from: https://github.com/kundajelab/tfmodisco/blob/178db48d9848d4d841ffda8c5d2f0269afa9f78e/modisco/visualization/viz_sequence.py
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def compute_per_position_ic(ppm, background, pseudocount):
    """
    Compute information content at each position of ppm.
    ppm: length x alphabet, rows sum to 1. background: base frequencies. pseudocount: added to ppm.
    Returns: 1d array of per-position information content.
    """
    assert len(ppm.shape) == 2
    assert ppm.shape[1] == len(background), "Letter axis must be second."
    if (not np.allclose(np.sum(ppm, axis=1), 1.0, atol=1.0e-5)):
        print("WARNING: Probabilities don't sum to 1 in all rows; renormalizing. PPM:\n" + str(ppm) + "\nSums:\n" + str(np.sum(ppm, axis=1)))
        ppm = ppm / np.sum(ppm, axis=1)[:, None]
    alphabet_len = len(background)
    ic = ((np.log((ppm + pseudocount) / (1 + pseudocount * alphabet_len)) / np.log(2))
          * ppm - (np.log(background) * background / np.log(2))[None, :])
    return np.sum(ic, axis=1)


def ic_scale(pwm, background):
    per_position_ic = compute_per_position_ic(ppm=pwm, background=background, pseudocount=0.001)
    return pwm * (per_position_ic[:, None])


def plot_a(ax, base, left_edge, height, color):
    a_polygon_coords = [
        np.array([[0.0, 0.0], [0.5, 1.0], [0.5, 0.8], [0.2, 0.0]]),
        np.array([[1.0, 0.0], [0.5, 1.0], [0.5, 0.8], [0.8, 0.0]]),
        np.array([[0.225, 0.45], [0.775, 0.45], [0.85, 0.3], [0.15, 0.3]])
    ]
    for polygon_coords in a_polygon_coords:
        ax.add_patch(matplotlib.patches.Polygon((np.array([1, height])[None, :] * polygon_coords + np.array([left_edge, base])[None, :]), facecolor=color, edgecolor=color))


def plot_c(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height, facecolor=color, edgecolor=color))
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height, facecolor='white', edgecolor='white'))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height, facecolor='white', edgecolor='white', fill=True))


def plot_g(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height, facecolor=color, edgecolor=color))
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height, facecolor='white', edgecolor='white'))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height, facecolor='white', edgecolor='white', fill=True))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.825, base + 0.085 * height], width=0.174, height=0.415 * height, facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.625, base + 0.35 * height], width=0.374, height=0.15 * height, facecolor=color, edgecolor=color, fill=True))


def plot_t(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.4, base], width=0.2, height=height, facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge, base + 0.8 * height], width=1.0, height=0.2 * height, facecolor=color, edgecolor=color, fill=True))


def plot_u(ax, base, left_edge, height, color):
    """
    Uracil (RNA channel 4). Do not use plot_t or thymine (T) is confused visually with U (ATGC illusion).
    Draw a U-shaped polygon in data coordinates like plot_a, width ~[0.18, 0.82] within the column,
    avoiding wide bold U from ax.text overlapping neighbors.
    """
    if abs(height) < 1e-12:
        return
    # Local coords x,y in [0,1]: outer outline + inner opening at top → U glyph (filled solid)
    xf = np.array(
        [0.18, 0.18, 0.28, 0.50, 0.72, 0.82, 0.82, 0.72, 0.72, 0.50, 0.28, 0.28],
        dtype=np.float64,
    )
    yf = np.array(
        [0.92, 0.28, 0.10, 0.04, 0.10, 0.28, 0.92, 0.92, 0.38, 0.18, 0.38, 0.92],
        dtype=np.float64,
    )
    poly = np.array([1, height])[None, :] * np.column_stack([xf, yf]) + np.array([left_edge, base])[None, :]
    ax.add_patch(
        matplotlib.patches.Polygon(
            poly, closed=True, facecolor=color, edgecolor=color, linewidth=0.0
        )
    )


default_colors = {0: 'green', 1: 'blue', 2: 'orange', 3: 'red'}
# RNA channel order A,C,G,U — channel 3 must be U, not plot_t
default_plot_funcs = {0: plot_a, 1: plot_c, 2: plot_g, 3: plot_u}


def plot_weights_given_ax(ax, array, figsize=(20, 2), height_padding_factor=0.2, length_padding=1.0, subticks_frequency=1.0, colors=default_colors, plot_funcs=default_plot_funcs, highlight={}, ylabel=""):
    if len(array.shape) == 3:
        array = np.squeeze(array)
    assert len(array.shape) == 2, array.shape
    if (array.shape[0] == 4 and array.shape[1] != 4):
        array = array.transpose(1, 0)
    assert array.shape[1] == 4
    array = np.asarray(array, dtype=np.float64)
    # Set y limits before drawing so plot_u transData heights match the final figure
    pos_stack = np.sum(np.maximum(array, 0.0), axis=1)
    neg_stack = np.sum(np.minimum(array, 0.0), axis=1)
    max_pos_height = float(np.max(pos_stack)) if pos_stack.size else 0.0
    min_neg_height = float(np.min(neg_stack)) if neg_stack.size else 0.0
    heights_at_positions = list(pos_stack)
    depths_at_positions = list(neg_stack)
    ax.set_xlim(-length_padding, array.shape[0] + length_padding)
    ax.xaxis.set_ticks(np.arange(0.0, array.shape[0] + 1, subticks_frequency))
    height_padding = max(abs(min_neg_height) * height_padding_factor, abs(max_pos_height) * height_padding_factor)
    ax.set_ylim(min_neg_height - height_padding, max_pos_height + height_padding)

    for i in range(array.shape[0]):
        per_base_vals = sorted(enumerate(array[i, :]), key=lambda x: abs(x[1]))
        positive_height_so_far = 0.0
        negative_height_so_far = 0.0
        for letter in per_base_vals:
            plot_func = plot_funcs[letter[0]]
            color = colors[letter[0]]
            if (letter[1] > 0):
                height_so_far = positive_height_so_far
                positive_height_so_far += letter[1]
            else:
                height_so_far = negative_height_so_far
                negative_height_so_far += letter[1]
            plot_func(ax=ax, base=height_so_far, left_edge=i, height=letter[1], color=color)
    for color in highlight:
        for start_pos, end_pos in highlight[color]:
            assert start_pos >= 0.0 and end_pos <= array.shape[0]
            min_depth = np.min(depths_at_positions[start_pos:end_pos])
            max_height = np.max(heights_at_positions[start_pos:end_pos])
            ax.add_patch(matplotlib.patches.Rectangle(xy=[start_pos, min_depth], width=end_pos - start_pos, height=max_height - min_depth, edgecolor=color, fill=False))
    ax.set_ylabel(ylabel)
    ax.yaxis.label.set_fontsize(15)


def plot_weights(array, figsize=(20, 2), save_path=None, **kwargs):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    plot_weights_given_ax(ax=ax, array=array, **kwargs)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()


def plot_weights_ax(array, ax, save_path=None, **kwargs):
    plot_weights_given_ax(ax=ax, array=array, **kwargs)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')


def plot_score_track_given_ax(arr, ax, threshold=None, **kwargs):
    ax.plot(np.arange(len(arr)), arr, **kwargs)
    if (threshold is not None):
        ax.plot([0, len(arr) - 1], [threshold, threshold])
    ax.set_xlim(0, len(arr) - 1)


def plot_score_track(arr, threshold=None, figsize=(20, 2), **kwargs):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    plot_score_track_given_ax(arr, threshold=threshold, ax=ax, **kwargs)
    plt.show()
