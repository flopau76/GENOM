import numpy as np
from collections import Counter
from abc import ABC, abstractmethod

from typing import Dict, List, Tuple

##########################################################################
#                              Metric class                              #
##########################################################################
""" This class only offers a slight packaging for different sliding windows. """

class Metric(ABC):
    @abstractmethod
    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        pass

    def update(self, **kwargs) -> None:
        for key, value in kwargs.items():
            setattr(self, key, value)


class distance(Metric):
    def __init__(self, ref_freq:Dict[int, float]=None, p:int=2) -> None:
        self.ref_freq = ref_freq
        self.p = p
        self.name = f"L{p}_distance"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        p = self.p
        ref_freq = self.ref_freq
        assert ref_freq is not None, "Reference frequency was not provided"

        kmers_count = Counter(kmers_list[:window_size])
        current_val = {kmer: abs(kmers_count[kmer]/window_size - ref_freq[kmer])**p for kmer in set(kmers_count.keys()).union(set(ref_freq.keys()))}

        current_res = sum(current_val.values())
        all_res = [current_res]

        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_res -= current_val[new_kmer] + current_val[old_kmer]
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            new_val = abs(kmers_count[new_kmer]/window_size - ref_freq[new_kmer])**p
            old_val = abs(kmers_count[old_kmer]/window_size - ref_freq[old_kmer])**p
            current_val[new_kmer] = new_val
            current_val[old_kmer] = old_val
            current_res += new_val + old_val

            if current_res < 0:
                current_res = 0
            all_res.append(current_res)
            
        return np.power(np.array(all_res),1/p)

class chi_squared(Metric):
    def __init__(self, ref_freq:Dict[int, float]=None, p:int=2) -> None:
        self.ref_freq = ref_freq
        self.name = "Chi-squared_distance"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        ref_freq = self.ref_freq
        assert ref_freq is not None, "Reference frequency was not provided"

        kmers_count = Counter(kmers_list[:window_size])
        current_val = {kmer: ((kmers_count[kmer]/window_size - ref_freq[kmer])**2)/ref_freq[kmer] for kmer in set(kmers_count.keys()).union(set(ref_freq.keys()))}

        current_res = sum(current_val.values())
        all_res = [current_res]

        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_res -= current_val[new_kmer] + current_val[old_kmer]
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            new_ref, old_ref = ref_freq[new_kmer], ref_freq[old_kmer]
            new_val = ((kmers_count[new_kmer]/window_size - new_ref)**2)/new_ref
            old_val = ((kmers_count[old_kmer]/window_size - old_ref)**2)/old_ref
            current_val[new_kmer], current_val[old_kmer] = new_val, old_val
            current_res += new_val + old_val
            all_res.append(current_res)
        return np.array(all_res)

class KLdivergence(Metric):
    def __init__(self, ref_freq:Dict[int, float]=None) -> None:
        self.ref_freq = ref_freq
        self.name = "KL_divergence"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        ref_freq = self.ref_freq
        assert ref_freq is not None, "Reference frequency was not provided"

        kmers_count = Counter(kmers_list[:window_size])
        current_freq = {kmer: count/window_size for kmer, count in kmers_count.items()}
        current_val = {kmer: p*np.log10(p/ref_freq[kmer]) for kmer, p in current_freq.items()}

        current_res = sum(current_val.values())
        all_res = [current_res]

        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_res -= current_val[old_kmer] + current_val.get(new_kmer,0)
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            new_freq = kmers_count[new_kmer]/window_size
            old_freq = kmers_count[old_kmer]/window_size
            new_val = new_freq*np.log10(new_freq/ref_freq[new_kmer])
            old_val = old_freq*np.log10(old_freq/ref_freq[old_kmer]) if old_freq > 0 else 0
            current_val[new_kmer] = new_val
            current_val[old_kmer] = old_val
            current_res += new_val + old_val
            all_res.append(current_res)
        return np.array(all_res)

class Jaccard(Metric):
    def __init__(self, ref_freq:Dict[int, int]=None) -> None:
        self.ref_freq = ref_freq
        self.name = "Jaccard_index"

    def slide_window(self, kmers_list:List[int], window_size:int) -> np.array:
        ref_freq = self.ref_freq
        assert ref_freq is not None, "Reference count was not provided"

        kmers_count = Counter(kmers_list[:window_size])
        current_inter = sum([min(count/window_size, ref_freq[kmer]) for kmer, count in kmers_count.items()])
        current_union = sum([max(count/window_size, ref_freq[kmer]) for kmer, count in kmers_count.items()])
        all_res = [current_inter/current_union]
        for i, new_kmer in  enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_inter -= min(kmers_count[old_kmer], ref_freq[old_kmer]) + min(kmers_count[new_kmer], ref_freq[new_kmer])
            current_union -= max(kmers_count[old_kmer], ref_freq[old_kmer]) + max(kmers_count[new_kmer], ref_freq[new_kmer])
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            current_inter += min(kmers_count[old_kmer], ref_freq[old_kmer]) + min(kmers_count[new_kmer], ref_freq[new_kmer])
            current_union += max(kmers_count[old_kmer], ref_freq[old_kmer]) + max(kmers_count[new_kmer], ref_freq[new_kmer])
            all_res.append(current_inter/current_union)
        return np.array(all_res)
    
class Convolution(Metric):
    def __init__(self, ref_count:Dict[int, int]=None) -> None:
        self.ref_count = ref_count
        self.name = "average_frequency"

    def slide_window(self, kmers_list:list[int], window_size : int) -> np.array:
        assert self.ref_count is not None, "Reference count was not provided"
        kmers_list = [self.ref_count[i] for i in kmers_list]
        all_avg = np.convolve(kmers_list, np.ones(window_size), 'valid')/window_size
        return -np.log10(all_avg)
