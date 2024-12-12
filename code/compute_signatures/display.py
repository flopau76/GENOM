from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf as bpdf

from typing import List

def display_windows(window_value:List[float], hits=None,
                     title:str=None, ylabel:str=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(window_value)
    if hits is not None:
        ax.scatter(hits.keys(), hits.values(), color='red', marker='+')
    ax.set_xlabel("Window start")
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    return fig

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


def save_pdf_fig_report(fig_list : List[plt.Axes], output_path_pdf):
    with bpdf.PdfPages(output_path_pdf) as pdf:
        for fig in fig_list:
            pdf.savefig(fig)
    return 0