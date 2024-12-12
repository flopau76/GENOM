import numpy as np
from collections import Counter
from abc import ABC, abstractmethod

from typing import Dict, List, Tuple

##########################################################################
#                               Metrics                                  #
##########################################################################

class Metric(ABC):
    """ This class is an abstract class for efficient metrics on sliding windows. 
    For it to work, we must have f(x1,...,xn) = h(g(x1)+...+g(xn)) 
    Having this, we can update f recomputing only the values g(xi) for which xi changed. """
    def __init__(self) -> None:
        self.name:str

    @abstractmethod
    def compute_single(self, kmer:int, count:int) -> float:
        """ Compute the metric for a single kmer. """
        pass

    def compute(self, kmers_count:Dict[int,int]):
        """ Compute the metric for a window. """
        return sum([self.compute_single(kmer, count) for kmer, count in kmers_count.items()])

    def slide_window(self, kmers_list:List[int], window_size:int) -> List[any]:
        """ Compute the metric for a whole genome. """
        kmers_count = Counter(kmers_list[:window_size])
        current_val = self.compute(kmers_count)
        all_val = [current_val]
        for i, new_kmer in enumerate(kmers_list[window_size:]):
            old_kmer = kmers_list[i]
            current_val = current_val - self.compute_single(new_kmer, kmers_count[new_kmer]) - self.compute_single(old_kmer, kmers_count[old_kmer])
            kmers_count[new_kmer] += 1
            kmers_count[old_kmer] -= 1
            current_val = current_val + self.compute_single(new_kmer, kmers_count[new_kmer]) + self.compute_single(old_kmer, kmers_count[old_kmer])
            all_val.append(current_val)
        return [self.compute_single(kmer, kmers_count[kmer]) for kmer in set(kmers_count.keys())]

    def post_process(self, value:List[any]) -> np.ndarray:
        return np.array(value)
    
    def compute_window(self, kmers_count:Dict[int,int]) -> float:
        return self.post_process(self.compute(kmers_count))[0]

##########################################################################
#                     Metric implementations                             #
##########################################################################

class total_value(Metric):
    def __init__(self, ref_value:Dict[int, float]) -> None:
        super().__init__()
        self.ref_value = ref_value
        self.name = "total_kmer_value"
    
    def compute_single(self, kmer:int, count:int) -> float:
        return count*self.ref_value[kmer]
    
##########################################################################

class distance(Metric):
    def __init__(self, ref_value:Dict[int, float], p:int=2, norm:float=1) -> None:
        super().__init__()
        self.ref_value = ref_value
        self.p = p
        self.norm = norm
        self.name = f"distance_{p}"

    def compute_single(self, kmer:int, count:int) -> float:
        return np.power(np.abs(count/self.norm - self.ref_value.get(kmer,0)), self.p)
    
    def compute(self, kmers_count:Dict[int,int]):
        return sum([self.compute_single(kmer, kmers_count[kmer]) for kmer in set(kmers_count.keys()).union(set(self.ref_value))])
    
    def post_process(self, value:float) -> float:
        return np.power(np.array(value), 1/self.p)
    
##########################################################################

class jaccard(Metric):
    def __init__(self, ref_value:Dict[int, float]) -> None:
        super().__init__()
        self.ref_value = ref_value
        self.name = "jaccard"

    def compute_single(self, kmer:int, count:int) -> float:
        inter = min(count, self.ref_value.get(kmer,0))
        union = max(count, self.ref_value.get(kmer,0))
        return np.array([inter, union])
    
    def compute(self, kmers_count:Counter[int]):
        inter = sum([min(count, self.ref_value.get(kmer,0)) for kmer, count in kmers_count.items()])
        union = sum([max(kmers_count[kmer], self.ref_value.get(kmer,0)) for kmer in set(kmers_count.keys()).union(set(self.ref_value.keys()))])
        return np.array([inter, union])

    def post_process(self, value:Tuple[int,int]) -> float:
        value = np.array(value)
        return value[:,0]/value[:,1]
    
##########################################################################

class KLdivergence(Metric):
    def __init__(self, ref_value:Dict[int, float], norm:float=1) -> None:
        super().__init__()
        self.ref_value = ref_value
        self.norm = norm
        self.name = "KLdivergence"

    def compute_single(self, kmer:int, count:int) -> float:
        def KLdiv(p, q):
            return p*np.log10(p/q) if p > 0 else 0
        return KLdiv(count/self.norm, self.ref_value.get(kmer,0))

##########################################################################
#       Old implementation using the kmer_list (for debugging only)      #
##########################################################################

def distance_kmers_list(kmers_list:List[int], ref_freq:Dict[int, float], window_size:int, p:int=2) -> np.array:
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

def KLdivergence_kmers_list(kmers_list:List[int], ref_freq:Dict[int, float], window_size:int) -> np.array:
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

def jaccard_kmers_list(kmers_list:List[int], ref_count:Dict[int, int], window_size:int) -> np.array:
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