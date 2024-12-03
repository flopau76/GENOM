from TP.loading import load_directory_as_pointers
from TP.kmers import stream_kmers_file
import TP.signatures as signatures
import numpy as np
import matplotlib.pyplot as plt
from time import time
from collections import Counter
import os, json

from typing import Dict, List


if __name__ == "__main__":
    k = 8
    input_folder = "data_test"
    output_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), f"transfer_summary_{input_folder}.json")
    
    nb_hits = 10
    best_hits = {}

    for sample, file_pointer in  load_directory_as_pointers(input_folder):
        start = time()
        print("Processing sample ", sample)
        print("    Computing the list of kmers")
        kmers_list = list(stream_kmers_file(file_pointer, k))

        print("    Starting frequence profile comparison")
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}

        window_average_values = signatures.window_slider_average(kmers_list, kmers_freq)
        window_distance = signatures.window_slider_euclidean_distance(kmers_list, kmers_freq, window_size=2000)

        # Get the nb_hits highest values of the window_distance array along with their indices
        # TODO: a pick spreads over several sliding windows. Instead, we should only look at local maxima
        highest_values_indices = np.argpartition(-window_distance, nb_hits)[:nb_hits]
        highest_values = window_distance[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}
        print("    Done in ", round(time()-start, 4), "s")

    #     fig, ax = plt.subplots()
    #     ax.plot(window_distance)
    #     ax.set_title("Sample "+sample)
    #     ax.set_xlabel("Window start index")
    #     ax.set_ylabel("Average kmer frequency")
    #     plt.show(block=False)
    #     plt.waitforbuttonpress(timeout=2)
    # plt.waitforbuttonpress()

    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)