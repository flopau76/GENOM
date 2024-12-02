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
        kmers_list = []
        dico = {}

        print("Processing", sample)
        size = len(file_pointer.read())
        file_pointer.seek(0)

        threshold = round(np.log2(size))
        
        kmers_list.append(list(stream_kmers_file(file_pointer, k)))

        kmer_split_list = np.array_split(*kmers_list, threshold)
        dico = dict(zip([f"{sample}_{i}" for i in range(threshold)], kmer_split_list))

        yield dico, kmers_list[0], sample


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


if __name__ == "__main__1":
    k = 8
    folder = "data_test"
    
    time_ = []
    out_dic = {}
    dir_path = '\\'.join(os.path.dirname(os.path.realpath(__file__)).split('\\')[:-1])

    
    for dico_strain, kmers_list, sample in compute_kmer(folder, k):
        st = time()
        liste_jac = compute_jaccard(dico_strain)
        df = km_stats.load_as_matrix(liste_jac)

        print("  Starting frequence profile comparison")
        dico_km = Counter(kmers_list)
        freq_dico_km = {key:n/sum(list(dico_km.values())) for key, n in Counter(kmers_list).items()}

        hits, freq_avg_km = km_stats.window_slider(kmers_list, freq_dico_km)
        time_.append(time()-st)

        out_dic[sample] = hits
        #disp.display_freq(freq_avg_km)

    with open(os.path.join(dir_path,"transfer_summary.json"), 'w') as outjson:
        json.dump(out_dic, outjson)
    
    print("Average runtime per genome", np.average(time_[1:]))
    print("Total Runtime", round(sum(time()-st)))
