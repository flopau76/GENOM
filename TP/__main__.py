from TP.loading import load_directory_as_pointers
from TP.kmers import stream_kmers_file
import numpy as np
from typing import Dict, List
from time import time
from collections import Counter
import os, json

import TP.km_stats as km_stats
import TP.display as disp

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
    folder = "data_test"
    
    time_ = []
    out_dic = {}

    dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    for sample, kmers_list in compute_kmer(folder, k):
        st = time()

        print("  Starting frequence profile comparison")
        dico_km = Counter(kmers_list)
        somme = sum(dico_km.values())
        freq_dico_km = {key:n/somme for key, n in dico_km.items()}

        hits, freq_avg_km = km_stats.window_slider(kmers_list, freq_dico_km)
        time_.append(time()-st)

        out_dic[sample] = hits
        disp.display_freq(freq_avg_km)
    
    with open(os.path.join(dir_path,"transfer_summary.json"), 'w') as outjson:
        json.dump(out_dic, outjson)
    
    print("Average runtime per genome", np.average(time_))
    print("Total Runtime", round(time()-st))
    