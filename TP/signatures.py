

import numpy as np
from collections import Counter
from scipy.signal import find_peaks
from scipy.special import kl_div

from typing import Dict

def window_slider_average(kmers_list:list[int], kmers_value:Dict[int, float], window_size:int = 2000) -> np.array:
    """Computes the average value of the kmers over each sliding"""
    current_val = np.sum([kmers_value[kmer] for kmer in kmers_list[0:window_size]])
    all_val = [current_val]

    for i in range(window_size, len(kmers_list)):
        current_val += kmers_value[kmers_list[i]] - kmers_value[kmers_list[i-window_size]]
        all_val.append(current_val)
    return np.array(all_val)/window_size

def window_slider_distance(kmers_list:list[int], kmers_ref_freq:Dict[int,float], p=2,  window_size:int=2000) -> np.array:
    """Computes the Lp distance between the signature of each sliding window and the reference signature
    Note 1: The signature of a sequence is the vector of kmer frequencies
    Note 2: The complexity of this function is O(nlog(n)) where n is the length of kmers_list"""
    current_count = Counter(kmers_list[0:window_size])
    diff_with_ref = {kmer:(current_count[kmer]/window_size - kmers_ref_freq[kmer]) for kmer in kmers_ref_freq}

    current_val = np.sum((diff)**p for diff in diff_with_ref.values())
    all_val = [current_val]
    for i in range(window_size, len(kmers_list)):
        current_val -= abs(diff_with_ref[kmers_list[i]])**p
        current_val -= abs(diff_with_ref[kmers_list[i-window_size]])**p
        current_count[kmers_list[i]] += 1
        current_count[kmers_list[i-window_size]] -= 1
        diff_with_ref[kmers_list[i]] = current_count[kmers_list[i]]/window_size - kmers_ref_freq[kmers_list[i]]
        diff_with_ref[kmers_list[i-window_size]] = current_count[kmers_list[i-window_size]]/window_size - kmers_ref_freq[kmers_list[i-window_size]]
        current_val += abs(diff_with_ref[kmers_list[i]])**p
        current_val += abs(diff_with_ref[kmers_list[i-window_size]])**p
        all_val.append(current_val)
    return np.power(np.array(all_val), 1/p)

def naive_KLdiv(kmers_list:list[int], kmers_ref_freq:Dict[int,float], window_size:int=2000):
    """
    Naive implementation of the Kullback-Leibler divergence to check
    output of the nlog(n) implementation.
    ==> Results are identical
    """
    all_div = []
    for i in range(0, len(kmers_list)-window_size+1):
        print(i)
        init_kmer_window = kmers_list[i:i+window_size]
        kl = 0
        for kmer in kmers_ref_freq.keys():
            freq_in_window = init_kmer_window.count(kmer)/window_size
            if freq_in_window > 0:
                kl += -(kmers_ref_freq[kmer]*np.log2(freq_in_window/kmers_ref_freq[kmer]))
        
        all_div.append(kl)
    return np.array(all_div)


def KLdivergence(kmers_list:list[int], kmers_ref_freq:Dict[int,float], window_size:int=2000):
    """
    Computes Kullback-Leibler divergence between two probability distribution.
    Here computes the divergence between the window k-mer probability distribution and the overall 
    k-mer distribution in the genome.
    """
    init_kmer_window = kmers_list[0:window_size]
    window_freq = {kmer: init_kmer_window.count(kmer)/window_size for kmer in kmers_ref_freq.keys()}

    current_kl_div = np.sum(-kmers_ref_freq[kmer]*np.log10(freq/kmers_ref_freq[kmer]) if freq > 0 else 0 for kmer, freq in window_freq.items())
    all_kl_div = [current_kl_div]

    for i in range(window_size, len(kmers_list)-window_size+1):
        first_kmer = kmers_list[i-window_size]
        P_first_kmer = kmers_ref_freq[kmers_list[i-window_size]]
        Q_first_kmer = window_freq[kmers_list[i-window_size]]

        new_kmer = kmers_list[i]
        P_new_kmer = kmers_ref_freq[new_kmer]
        Q_new_kmer = window_freq[new_kmer]

        current_kl_div -= -(P_first_kmer*np.log2(Q_first_kmer/P_first_kmer))
        if Q_new_kmer > 0: #avoid warning for log(0)
            current_kl_div -= -(P_new_kmer*np.log2(Q_new_kmer/P_new_kmer))

        window_freq[first_kmer] -= 1/window_size
        window_freq[new_kmer] += 1/window_size

        if window_freq[first_kmer] > 0: #avoid warning for log(0)
            current_kl_div += -(P_first_kmer*np.log2(window_freq[first_kmer]/P_first_kmer))
        current_kl_div += -(P_new_kmer*np.log2(window_freq[new_kmer]/P_new_kmer))

        all_kl_div.append(current_kl_div)   
    return np.array(all_kl_div)

def find_maxima(window_array:np.array, nb_maxima:int) -> np.array:
    """Finds the indices of the nb_maxima highest values in the window_array"""
    return np.argpartition(-window_array, nb_maxima)[:nb_maxima]

def find_local_maxima(window_array:np.array, **kwargs) -> np.array:
    """Finds the indices of the local maxima above a given threshold in the window_array"""
    return find_peaks(window_array, **kwargs)[0]