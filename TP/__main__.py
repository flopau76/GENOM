import numpy as np
import matplotlib.pyplot as plt

from collections import Counter
from itertools import product
import matplotlib.pyplot as plt

import TP.window_slider as window_slider
from TP.loading import iter_directory, open_genome
import TP.display as display

from scipy.signal import find_peaks
from time import time
import os, json

from typing import Dict, List

def dump_matrix_to_csv(matrix, filename, delimiter=',', precision=6):
    """ Dump a NumPy matrix to a CSV file. """
    
    if not filename.lower().endswith('.csv'):
        filename += '.csv'
    
    try:
        np.set_printoptions(precision=precision, suppress=True)
        np.savetxt(filename, matrix, delimiter=delimiter, fmt=f'%0.{precision}f')
        print(f"Matrix successfully exported to: {filename}")
    
    except Exception as e:
        print(f"Error exporting matrix to CSV: {e}")

if __name__ == "__main__":
    k = 8
    window_size = 100
    step = 5
    input_folder = "toy_transfer"
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    files = iter_directory(input_folder)
    sample_a, file_a = files[0]
    sample_b, file_b = files[1]

    file_a = open_genome(file_a)
    file_b = open_genome(file_b)

    res = []
    start = time()
    for window_b in window_slider.stream_windows(file_b, window_size, k, step):
        res.append(window_slider.compute_metrics_file(file_a, ["jacc"], window_size, k, step=step, ref_count=window_b.count))
    res = np.array(res)
    print(f"Done in: {time() - start:.2f}s")

    display.display_matrix(res, f"test_ws{window_size}_s{step}.png")
    plt.show(block=True)
