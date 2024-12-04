from TP.loading import load_directory_as_pointers
from TP.kmers import stream_kmers_file
import TP.signatures as signatures

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.signal import find_peaks
from time import time
import os, json

from typing import Dict, List


if __name__ == "__main__":
    k = 8
    window_size = 2000
    input_folder = "data_test"
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    best_hits = {}

    for sample, file_pointer in  load_directory_as_pointers(input_folder):
        start = time()
        print("Processing sample ", sample)
        print("    Computing the list of kmers")
        kmers_list = list(stream_kmers_file(file_pointer, k))

        print("    Starting frequence profile comparison")
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}
        window_distance = signatures.window_slider_distance(kmers_list, kmers_freq, window_size=window_size)

        # kmers_rarity = {kmer:1/freq for kmer, freq in kmers_freq.items()}
        # window_average_rarity = signatures.window_slider_average(kmers_list, kmers_rarity, window_size)

        print("    Finding the best hits")
        # highest_values_indices = signatures.find_maxima(window_distance, nb_hits)
        # highest_values = window_distance[highest_values_indices]

        height = np.mean(window_distance) + 10*np.var(window_distance)
        distance = window_size
        prominence = (np.max(window_distance) - np.min(window_distance))/3
        highest_values_indices, _ = find_peaks(window_distance, prominence=prominence)     # note: dependig on the metric, filtering parameters may need to be adjusted
        highest_values = window_distance[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

        print("    Done in ", round(time()-start, 4), "s")

    #     # Plotting
    #     fig, ax = plt.subplots()
    #     ax.plot(window_distance)
    #     ax.plot(highest_values_indices, highest_values, 'r*')
    #     ax.set_title("Sample "+sample)
    #     ax.set_xlabel("Window start index")
    #     ax.set_ylabel("Distance between window and average signature")
    #     plt.show(block=False)
    #     plt.waitforbuttonpress(timeout=2)
    # plt.waitforbuttonpress()

    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)