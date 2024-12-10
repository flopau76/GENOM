from TP.loading import load_directory_as_pointers
from TP.kmers import stream_kmers_file

import TP.signatures as signatures
import TP.km_stats as km_stats
import TP.display as disp

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.signal import find_peaks
from time import time
import os, json
from typing import Dict, List

from typing import Dict, List

def dict_intersection(dictA, dictB):
    """ Computes the intersection of two dictionaries
    :param dict dictA, dictB: dictionaries to compare"""
    intersection = 0
    for key, val in dictA:
        intersection += min(val, dictB.get(key, 0))
    return intersection

def list_intersection(listA, listB):
    """ Computes the intersection of two sorted lists
    :param np.array listA, listB: sorted np.array to compare"""
    intersection = 0
    idxA = 0
    idxB = 0
    while idxA < len(listA) and idxB < len(listB):
        if listA[idxA] == listB[idxB]:
            intersection += 1
            idxA += 1
            idxB += 1
        elif listA[idxA] < listB[idxB]:
            idxA += 1
        else:
            idxB += 1
    return intersection

def xorshift(val):
    """ Hash function using the xorshift algorithm """
    val ^= val << 13
    val &= 0xFFFFFFFFFFFFFFFF
    val ^= val >> 7
    val ^= val << 17
    val &= 0xFFFFFFFFFFFFFFFF
    return val

def compute_kmer(folder : str, k : int):
    # Computing the kmers
    print("  Computing the kmers")
    for sample, file_pointer in load_directory_as_pointers(folder):
        print("Processing", sample)        
        yield sample, list(stream_kmers_file(file_pointer, k))

def compute_jaccard(dico : Dict[str, List[int]]):
    # Computing the Jaccard index
    print("  Computing the pairwise similarities")
    filenames = list(dico.keys())
    list_tuple_jac = []

    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            intersection = list_intersection(sorted(dico[filenames[i]].tolist()), sorted(dico[filenames[j]].tolist()))
            dist_j = intersection / (len(dico[filenames[i]]) + len(dico[filenames[j]]) - intersection)
            #print(f"{'==='*20}\n{filenames[i]} | {filenames[j]} | {dist_j}")

            list_tuple_jac.append((filenames[i], filenames[j], dist_j))

    return list_tuple_jac


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
        
        st = time()
        window_distance = signatures.window_slider_distance(kmers_list, kmers_freq, window_size=window_size)
        #disp.display_windows(window_distance, ylabel="L2 distance", title="L2 distance for sample "+sample)
        print("   L2 done in ", round(time()-st, 4), "s")

        st = time()
        window_kldiv = signatures.KLdivergence(kmers_list, kmers_freq)
        print("   KL divergencce done in ", round(time()-st, 4), "s")
        #disp.display_windows(window_kldiv, ylabel="KL divergence", title="KL divergence for sample "+sample)
        
        #kldiv2 = signatures.naive_KLdiv(kmers_list, kmers_freq)
        #disp.display_freq(kldiv2)
        
        # kmers_rarity = {kmer:1/freq for kmer, freq in kmers_freq.items()}
        # window_average_rarity = signatures.window_slider_average(kmers_list, kmers_rarity, window_size)

        
        #prop_test = window_kldiv/window_average_rarity
        #plt.plot(list(range(len(window_average_rarity))), prop_test)        
        
        print("    Finding the best hits")
        # highest_values_indices = signatures.find_maxima(window_distance, nb_hits)
        # highest_values = window_distance[highest_values_indices]


        """
        # possible filtering parameters: dependig on the metric, they may need to be adjusted
        height = np.mean(window_distance) + 10*np.var(window_distance)
        distance = window_size"""
        prominence = (np.max(window_distance) - np.min(window_distance))/3
        highest_values_indices, _ = find_peaks(window_distance, prominence=prominence)
        highest_values = window_distance[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

        print("    Done in ", round(time()-start, 4), "s")

    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)
