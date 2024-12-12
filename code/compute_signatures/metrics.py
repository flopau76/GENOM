import numpy as np
from collections import Counter
from abc import ABC, abstractmethod

from typing import Dict, List, Tuple

##########################################################################
#       Old implementation using the kmer_list (for debugging only)      #
##########################################################################

class Metric(ABC):
    @abstractmethod
    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        pass

class distance(Metric):
    def __init__(self, ref_value:Dict[int, float], p:int=2, norm:float=1) -> None:
        self.ref_value = ref_value
        self.p = p
        self.name = f"L{p} distance"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        p = self.p
        ref_freq = self.ref_value

        kmers_count = Counter(kmers_list[:window_size])
        current_diff = {kmer: kmers_count[kmer]/window_size - ref_freq[kmer] for kmer in set(kmers_count.keys()).union(set(ref_freq.keys()))}
        current_val = sum([abs(diff)**p for diff in current_diff.values()])
        all_val = [current_val]

        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_val -= abs(current_diff[new_kmer])**p
            current_val -= abs(current_diff[old_kmer])**p
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            current_diff[new_kmer] = kmers_count[new_kmer]/window_size - ref_freq[new_kmer]
            current_diff[old_kmer] = kmers_count[old_kmer]/window_size - ref_freq[old_kmer]
            current_val += abs(current_diff[new_kmer])**p
            current_val += abs(current_diff[old_kmer])**p
            all_val.append(current_val)
        return np.power(np.array(all_val),1/p)
    
class KLdivergence(Metric):
    def __init__(self, ref_freq:Dict[int, float]) -> None:
        self.ref_freq = ref_freq
        self.name = "KL divergence"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        ref_freq = self.ref_freq

        kmers_count = Counter(kmers_list[:window_size])
        current_freq = {kmer: count/window_size for kmer, count in kmers_count.items()}
        current_val = sum([p*np.log10(p/ref_freq[kmer]) for kmer, p in current_freq.items()])
        all_val = [current_val]

        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            new_p = current_freq[old_kmer]
            old_p = current_freq[new_kmer]
            current_val -= new_p*np.log10(new_p/ref_freq[new_kmer])
            current_val -= old_p*np.log10(old_p/ref_freq[old_kmer])
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            current_freq[new_kmer] = kmers_count[new_kmer]/window_size
            current_freq[old_kmer] = kmers_count[old_kmer]/window_size
            new_p = current_freq[old_kmer]
            old_p = current_freq[new_kmer]
            current_val += old_p*np.log10(old_p/ref_freq[new_kmer]) if old_p > 0 else 0
            current_val += new_p*np.log10(new_p/ref_freq[old_kmer])
        return np.array(all_val)

class Jaccard(Metric):
    def __init__(self, ref_count:Dict[int, int]) -> None:
        self.ref_count = ref_count
        self.name = "Jaccard index"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        ref_count = self.ref_count

        kmers_count = Counter(kmers_list[:window_size])
        current_inter = sum([min(count, ref_count[kmer]) for kmer, count in kmers_count.items()])
        current_union = sum([max(count, ref_count[kmer]) for kmer, count in kmers_count.items()])
        all_val = [current_inter/current_union]
        for i, new_kmer in  enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_inter -= min(kmers_count[old_kmer], ref_count[old_kmer]) + min(kmers_count[new_kmer], ref_count[new_kmer])
            current_union -= max(kmers_count[old_kmer], ref_count[old_kmer]) + max(kmers_count[new_kmer], ref_count[new_kmer])
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            current_inter += min(kmers_count[old_kmer], ref_count[old_kmer]) + min(kmers_count[new_kmer], ref_count[new_kmer])
            current_union += max(kmers_count[old_kmer], ref_count[old_kmer]) + max(kmers_count[new_kmer], ref_count[new_kmer])
            all_val.append(current_inter/current_union)
        return np.array(all_val)