import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from collections import Counter

from time import time
import os, json

import compute_signatures.metrics as metrics
from compute_signatures.loading import iter_directory, open_genome
from compute_signatures.kmers import stream_kmers_file
from compute_signatures.display import display_windows, save_pdf_fig_report

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """
    Show a progress bar indicating downloading progress
    """
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()


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
    folder_name = "toy_transfer"

    base_dir =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    input_folder = os.path.join(base_dir, "input", "sequence_db", folder_name)
    output_path = os.path.join(base_dir, "output", f"transfer_summary_{folder_name}.json")
    output_path_pdf = os.path.join(base_dir, "output", f"transfer_summary_{folder_name}.pdf")
    best_hits = {}

    fig_list = []
    n=1

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
        fig_list.append(fig)

        # finding best hits:
        distances = results[0,:]
        highest_values_indices, _ = find_peaks(distances, prominence=(np.max(distances) - np.min(distances))/3)
        highest_values = distances[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

        progressbar(n, len(os.listdir(input_folder)))
        n+=1

    save_pdf_fig_report(fig_list, output_path_pdf)
    # saving hits to json
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)

    # waiting for all figures to be closed manually
    #plt.show(block=True)