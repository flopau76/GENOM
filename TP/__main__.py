from TP.loading import load_directory, load_directory_as_pointers
from TP.kmers import stream_kmers, stream_kmers_file
import numpy as np
from typing import Dict, List

from TP.km_stats import load_as_matrix

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

        threshold = 50#round(np.log2(size))
        kmers_list.append(list(stream_kmers_file(file_pointer, k)))

        kmer_split_list = np.array_split(*kmers_list, threshold)
        dico = dict(zip([f"{sample}_{i}" for i in range(threshold)], kmer_split_list))

        yield dico, threshold

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

# Faire rapport statistic global 

if __name__ == "__main__":
    k = 8
    folder = "data_test"

    for dico_strain, threshold in compute_kmer(folder, k):
        liste_jac = compute_jaccard(dico_strain)
        print(load_as_matrix(liste_jac, threshold))
        break
