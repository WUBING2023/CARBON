"""Loss curve and basic training plots."""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


def plotLoss(y_loss, plotName, bin=False, binSize=1):
    if bin:
        y_loss = [np.mean(y_loss[i:i + binSize]) for i in range(0, len(y_loss), binSize)]
    x_range = range(len(y_loss))
    plt.figure()
    ax = plt.axes()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.plot(x_range, y_loss, color='blue', linewidth=1.5, linestyle="solid")
    plt.legend(['Train Loss'], loc='upper right')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Number of training examples')
    plt.ylabel('Loss')
    plt.title('Loss curve')
    plt.xlim((0, len(y_loss)))
    plt.savefig(f"{plotName}_{binSize}_loss.png")
    plt.show()
