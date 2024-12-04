

import numpy as np
from collections import Counter
from scipy.signal import find_peaks

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

def find_maxima(window_array:np.array, nb_maxima:int) -> np.array:
    """Finds the indices of the nb_maxima highest values in the window_array"""
    return np.argpartition(-window_array, nb_maxima)[:nb_maxima]

def find_local_maxima(window_array:np.array, **kwargs) -> np.array:
    """Finds the indices of the local maxima above a given threshold in the window_array"""
    return find_peaks(window_array, **kwargs)[0]