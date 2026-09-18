"""Attribution logo and figure plotting (logomaker)."""

import logomaker
import pandas as pd
import numpy as np


def plot_sequence_attribution(attribution_matrix, bases=['A', 'C', 'G', 'T'], width=20, height=4):
    attribution_df = pd.DataFrame(attribution_matrix, columns=bases)
    attribution_logo = logomaker.Logo(attribution_df, shade_below=.5, fade_below=.5, font_name='Arial Rounded MT Bold')
    attribution_logo.style_spines(visible=False)
    attribution_logo.style_spines(spines=['left', 'bottom'], visible=True)
    attribution_logo.ax.set_ylabel("IG Attribution", labelpad=-1)
    attribution_logo.fig.set_figheight(height)
    attribution_logo.fig.set_figwidth(width)


def make_attribution_figure(contribution_score, ax):
    df = pd.DataFrame(contribution_score, columns=['A', 'C', 'G', 'U'])
    logo = logomaker.Logo(df, shade_below=0, fade_below=0, font_name='Arial Rounded MT Bold', ax=ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_xticks(np.arange(0, df.shape[0], 100))
    ax.set_xlim(0, df.shape[0])
    ax.set_yticks([])
    ax.set_ylim(np.min(contribution_score) - 0.25, 1)
    ax.relim()
    ax.autoscale_view()


def make_attribution_logomaker(contribution_score, ax):
    df = pd.DataFrame(contribution_score, columns=['A', 'C', 'G', 'U'])
    logo = logomaker.Logo(df, shade_below=0, fade_below=0, ax=ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks([])
    ax.set_xticks(np.arange(0, df.shape[0], 100))
    ax.relim()
    ax.autoscale_view()
