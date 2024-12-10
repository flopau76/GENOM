import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from collections import Counter

from time import time
import os, json

import TP.metrics as metrics
from TP.loading import iter_directory, open_genome
from TP.kmers import stream_kmers_file
from TP.display import display_windows

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
    window_size = 2000
    input_folder = "toy_transfer"
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    best_hits = {}

    for sample, file_name in  iter_directory(input_folder):
        file_pointer = open_genome(file_name)

        # iterate once over the file to compute the total kmer frequency
        kmers_list = list(stream_kmers_file(file_pointer, k))
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}

        # compare each window to global kmer frequency using different metrics
        metric_list = [metrics.distance(kmers_freq, norm=window_size), metrics.KLdivergence(kmers_freq, norm=window_size)]
        results = metrics.compute_metrics_file(file_pointer, metric_list, window_size, k)
        file_pointer.close()

        # display results
        fig = plt.figure()
        for i, metric in enumerate(metric_list):
            ax = fig.add_subplot(len(metric_list), 1, i+1)
            ax.plot(results[i,:])
            ax.set_title(metric.name)
        ax.set_xlabel("Position")
        fig.suptitle(f"Metrics for sample {sample}")
        fig.tight_layout()
        # fig.savefig(f"metrics_{sample}.png")

        # finding best hits:
        distances = results[0,:]
        highest_values_indices, _ = find_peaks(distances, prominence=(np.max(distances) - np.min(distances))/3)
        highest_values = distances[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

    # saving hits to json
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)

    # waiting for all figures to be closed manually
    plt.show(block=True)