from matplotlib import pyplot as plt

from typing import List

def display_windows(window_value:List[float], title:str=None, ylabel:str=None, ax=None, label:str=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(window_value, label=label)
    ax.set_xlabel("Window position")
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if label is not None:
        ax.legend()
    plt.show(block=False)
    plt.waitforbuttonpress(timeout=0.5)
    return 0

def display_matrix(matrix, save_path, step=1, title='Matrix Heatmap', 
                     figsize=(10,10), cmap='viridis', colorbar=True):
    """ Visualize a NumPy matrix as a heatmap using Matplotlib. """
    fig, ax = plt.subplots(figsize=figsize)
    extent = [0, matrix.shape[1]*step, 0, matrix.shape[0]*step]
    im = ax.imshow(matrix, cmap=cmap, aspect='equal',interpolation ='none', extent=extent)
    ax.set_title(title)

    if colorbar:
        plt.colorbar(im, ax=ax, label='Value')
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=1000)
    return fig