import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from collections import Counter
import argparse, sys
import os, json

import compute_signatures.metrics as metrics
from compute_signatures.loading import iter_directory, open_genome
from compute_signatures.kmers import stream_kmers_file

from typing import List, Dict
from io import TextIOWrapper

from time import time

def find_potential_HGT(file_pointer:TextIOWrapper, window_size:int, k:int, metric:metrics.Metric, step:int=1) -> Dict[int, float]:
    """ Find potential HGT regions in a genome file """
    # slide over the file to compute the metric
    results = metrics.compute_metrics_file(file_pointer, metric_list, window_size, k, step)[0]

    # find best hits:
    highest_values_indices, _ = find_peaks(results, prominence=(np.max(results) - np.min(results))/3)
    highest_values = results[highest_values_indices]
    best_hits = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

    return best_hits

def plot_profile(file_pointer:TextIOWrapper, window_size:int, k:int, metric_list:List[metrics.Metric], step:int=1):
    """ Plot the value of some metrics along a genome """
    # slide over the file to compute the metric
    results = metrics.compute_metrics_file(file_pointer, metric_list, window_size, k, step)[0]

    # display results
    fig = plt.figure()
    for i, metric in enumerate(metric_list):
        ax = fig.add_subplot(len(metric_list), 1, i+1)
        ax.plot(results[i])
        ax.set_title(metric.name)
    ax.set_xlabel("Position")
    fig.suptitle(f"Metrics for sample {sample}")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    k = 8
    window_size = 2000
    input_folder = "procariote_diversity"
    # input_folder = "toy_transfer"

    input_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../../input/sequence_db", input_folder)
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    best_hits = {}

    for sample, file_name in  iter_directory(input_path):
        file_pointer = open_genome(file_name)

        # iterate once over the file to get the global signature
        t = time()
        kmers_list = list(stream_kmers_file(file_pointer, k))
        print(f"Time for computing the kmers: {time()-t}")
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}
        print(f"Sample {sample} has {len(kmers_freq)} kmers")

        # compute metrics using the metric class
        t0 = time()
        metric_list = [metrics.distance(kmers_freq, norm=window_size)]
        kmers_stream = stream_kmers_file(file_pointer, k)
        result = metrics.compute_metrics_file(kmers_stream, metric_list, window_size)[0]
        print(f"Time for distance: {time()-t0}")

        # compute metrics using the old function
        t1 = time()
        result_1 = metrics.distance_kmers_list(kmers_list, kmers_freq, window_size, p=2)
        print(f"Time for distance_kmers_list: {time()-t1}")

        t2 = time()
        results_2 = metrics.window_slider_distance(kmers_list, kmers_freq, window_size, p=2)
        print(f"Time for window_slider_distance: {time()-t2}")

        plt.plot(result, label="new")
        plt.plot(result_1, label="old Flo")
        plt.plot(results_2, label="old Matt")
        plt.legend()
        plt.show()

        print("Done computing distance")

    # waiting for all figures to be closed manually
    plt.show(block=True)




if __name__ == "__main__new":

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--display', help='Display the plots', action='store_true')
    parser.add_argument("input_db", help="The name of the input database (must be in `input/sequence_db/`)")
    parser.add_argument('-k', '--kmer', help='The size of the kmer', type=int, default=8)
    parser.add_argument('-w', '--window', help='The size of the window', type=int, default=2000)

    args = parser.parse_args()

    root_dir = os.path.abspath(os.path.join(__file__, "..", "..", ".."))
    input_folder = os.path.join(root_dir, "input", "sequence_db", args.input_db)
    output_path = os.path.join(root_dir, "output" f"transfer_summary_{input_folder}.json")

    best_hits = dict()

    for sample, file_name in  iter_directory(input_folder):
        file_pointer = open_genome(file_name)

        # iterate once over the file to compute the total kmer frequency
        kmers_count = Counter(stream_kmers_file(file_pointer, k))
        kmers_tot = sum(kmers_count.values())
        kmers_freq = {kmer:count/kmers_tot for kmer, count in kmers_count.items()}

        # use the global signature to create some metrics
        metric_list = [metrics.distance(kmers_freq, norm=window_size), metrics.KLdivergence(kmers_freq, norm=window_size)]

        # find potential HGT regions
        best_hits[sample] = find_potential_HGT(file_pointer, window_size, k, metric_list[0])

        # plot the profile of the metrics
        fig = plot_profile(file_pointer, window_size, k, metric_list)

    # saving hits to json
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)

    # waiting for all figures to be closed manually
    plt.show(block=True)