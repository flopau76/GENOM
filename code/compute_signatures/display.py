from matplotlib import pyplot as plt

from collections import defaultdict
from matplotlib.patches import Rectangle

from typing import List
import numpy as np
import os, json

def parse_HGT_report(file:str):
    res = {}
    with open(file, 'r') as f:
        f.readline()
        for line in f.readlines():
            line = line.split('\t\t')
            send_start, send_end, rec_name, rec_start = line[1:]
            send_start, send_end, rec_start = int(send_start), int(send_end), int(rec_start)
            res[rec_name]= {'pos':[[rec_start], [rec_start+send_end-send_start]]}
    return res

def get_ground_truth(folder:str):
    """" Looks for a file called ground_truth.json or HGT_report.txt in the folder
    and transforms it into a dictionnary containing the ground truth """

    HGT_path = os.path.join(folder, "HGT_report.txt")
    if os.path.exists(HGT_path):
        return parse_HGT_report(HGT_path)
    
    ground_truth_path = os.path.join(folder, "ground_truth.json")
    if os.path.exists(ground_truth_path):
        with open(ground_truth_path, 'r') as gt_file:
            ground_truth = json.load(gt_file)
        return ground_truth
    else:
        return defaultdict(None)

def display_windows(window_value:List[float], hits=None, ground_truth=None,
                     title:str=None, ylabel:str=None, dpi=None) -> plt.Figure:
    fig, ax = plt.subplots(dpi=dpi)
    ax.plot(window_value)
    if hits is not None:
        ax.scatter(hits.keys(), hits.values(), color='red', marker='+')
    if ground_truth is not None:
        ground_truth_pos, ground_truth_neg = ground_truth.get('pos', None), ground_truth.get('neg', None)
        if ground_truth_pos is not None:
            rect = Rectangle((ground_truth_pos[1][0], min(window_value)), width=(ground_truth_pos[0][0]-ground_truth_pos[1][0])*2, height=max(window_value)-min(window_value), edgecolor='red', facecolor="white", alpha = 1, label='Known HGT')
            ax.add_patch(rect)

        if ground_truth_neg is not None:
            rect = Rectangle((ground_truth_neg[1][0], min(window_value)), width=(ground_truth_pos[0][0]-ground_truth_pos[1][0])*2, height=max(window_value)-min(window_value), edgecolor='red', facecolor="white", alpha = 1, label='Known HGT')
            ax.add_patch(rect)

    ax.set_xlabel("Window start")
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if ground_truth is not None:
        plt.legend()

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