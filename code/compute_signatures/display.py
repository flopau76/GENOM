##########################################################################
#                        Useless/Deprecated ?                            #
##########################################################################


from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
from matplotlib.patches import Rectangle

from typing import List, Dict
import numpy as np

def ref_parse(ref_path : str):
    with open(ref_path, 'r') as ref_file:
        ref_file =[line.split('\t') for line in ref_file.readlines()]

    res = {line[4]:int(line[6].strip()) for index, line in enumerate(ref_file) if index != 0}
    return res    

def display_windows(window_value:List[float], sample : str, hits=None,
                     title:str=None, ylabel:str=None, ax=None, 
                     dpi=None, ref : Dict[str,int] = None, window_size : int = 5000) -> plt.Figure: 

    if ax is None:
        fig, ax = plt.subplots(dpi=dpi)
    ax.plot(window_value)
    if ref is not None:
        try :
            known_HGT_position = ref[sample]
            size = (len(window_value)/window_size)*np.log2(len(window_value))*np.log10(window_size)
            rect = Rectangle((known_HGT_position-(size/2), min(window_value)), width=size, height=max(window_value)-min(window_value), edgecolor='red', facecolor="white", label="Known HGT")
            ax.add_patch(rect)
            ax.legend(loc='center left', bbox_to_anchor=(0.8, 1))
        except Exception:
            pass
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