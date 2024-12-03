

import numpy as np
from collections import Counter

from typing import Dict

def window_slider_average(kmers_list:list[int], kmers_value:Dict[int, float], window_size:int = 2000) -> np.array:
    """Computes the average value of the kmers over each sliding"""
    current_val = np.sum([kmers_value[kmer] for kmer in kmers_list[0:window_size]])
    all_val = [current_val]

    for i in range(window_size, len(kmers_list)):
        current_val += kmers_value[kmers_list[i]] - kmers_value[kmers_list[i-window_size]]
        all_val.append(current_val)
    return np.array(all_val)/window_size

def window_slider_euclidean_distance(kmers_list:list[int], kmers_ref_freq:Dict[int,float],  window_size:int=2000) -> np.array:
    """Computes the euclidean distance between the signature of each sliding window and the reference signature
    Note 1: The signature of a sequence is the vector of kmer frequencies
    Note 2: The complexity of this function is O(nlog(n)) where n is the length of kmers_list"""
    current_count = Counter(kmers_list[0:window_size])
    diff_with_ref = {kmer:(current_count[kmer]/window_size - kmers_ref_freq[kmer]) for kmer in kmers_ref_freq}

    current_val = np.sum((diff)**2 for diff in diff_with_ref.values())
    all_val = [current_val]
    for i in range(window_size, len(kmers_list)):
        current_val -= (diff_with_ref[kmers_list[i]])**2
        current_val -= (diff_with_ref[kmers_list[i-window_size]])**2
        current_count[kmers_list[i]] += 1
        current_count[kmers_list[i-window_size]] -= 1
        diff_with_ref[kmers_list[i]] = current_count[kmers_list[i]]/window_size - kmers_ref_freq[kmers_list[i]]
        diff_with_ref[kmers_list[i-window_size]] = current_count[kmers_list[i-window_size]]/window_size - kmers_ref_freq[kmers_list[i-window_size]]
        current_val += (diff_with_ref[kmers_list[i]])**2
        current_val += (diff_with_ref[kmers_list[i-window_size]])**2
        all_val.append(current_val)
    return np.power(np.array(all_val), 0.5)