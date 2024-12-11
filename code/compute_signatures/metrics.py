import numpy as np
from collections import Counter
from abc import ABC, abstractmethod
from itertools import islice

from compute_signatures.window_slider import window

from typing import Dict, List, Tuple, Generator

##########################################################################
#                               Metrics                                  #
##########################################################################

def get_subset(d:Dict[int,int], changes:Dict[int, int]) -> Dict[int,int]:
    """ Returns the subset of d that is affected by the changes """
    old = dict()
    current = dict()
    for k, c in changes.items():
        old[k] = d[k] - c
        current[k] = d[k]
    return old, current

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

    def update(self, old_value, new_count:Dict[int,int], new_kmer, old_kmer) -> float:
        """ Update the metric value after a window slide """
        count_new = new_count[new_kmer]
        count_old = new_count[old_kmer]
        return old_value + self.compute_single(new_kmer, count_new) + self.compute_single(old_kmer, count_old) \
                        - self.compute_single(new_kmer, count_new-1) - self.compute_single(old_kmer, count_old+1)

    def post_process(self, value:List[any]) -> np.ndarray:
        return np.array(value)
    
    def compute_window(self, kmers_count:Dict[int,int]) -> float:
        return self.post_process(self.compute(kmers_count))[0]

##########################################################################
#                Iterate over file to compute metric                     #
##########################################################################

def compute_metrics_file(kmers_stream, metric_list:List[Metric], window_size:int) -> List[np.ndarray]:
    """ Compute the metric for all windows of a file
    Currently supported metrics are: "average", "distance", "KLdiv", "jacc" """

    current_window = window(islice(kmers_stream, window_size))
    current_results = [metric.compute(current_window.count) for metric in metric_list]
    results = [current_results]

    for new_kmer in kmers_stream:
        old_kmer = current_window.slide_right(new_kmer)
        current_results = [metric.update(old_value, current_window.count, new_kmer, old_kmer) for metric, old_value in zip(metric_list, current_results)]
        results.append(current_results)

    new_results = []
    for i, metric in enumerate(metric_list):
        res = [result[i] for result in results]
        new_results.append(metric.post_process(res))
    return new_results

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
    current_val = sum([min(count, ref_count[kmer])/max(count, ref_count[kmer]) for kmer, count in kmers_count.items()])
    all_val = [current_val]
    for i, new_kmer in  enumerate(kmers_list[window_size:]):
        old_kmer = kmers_list[i]
        current_val -= min(kmers_count[old_kmer], ref_count[old_kmer])/max(kmers_count[old_kmer], ref_count[old_kmer])
        current_val -= min(kmers_count[new_kmer], ref_count[new_kmer])/max(kmers_count[new_kmer], ref_count[new_kmer])
        kmers_count[new_kmer] += 1
        kmers_count[old_kmer] -= 1
        current_val += min(kmers_count[new_kmer], ref_count[new_kmer])/max(kmers_count[new_kmer], ref_count[new_kmer])
        current_val += min(kmers_count[old_kmer], ref_count[old_kmer])/max(kmers_count[old_kmer], ref_count[old_kmer])
        all_val.append(current_val)
    return np.array(all_val)


##########################################################################
#            Old implementation 2 (for debugging only)                   #
##########################################################################


def window_slider_average(kmers_list:list[int], kmers_value:Dict[int, float], window_size:int = 2000) -> np.array:
    """Computes the average value of the kmers over each sliding"""
    current_val = np.sum([kmers_value[kmer] for kmer in kmers_list[0:window_size]])
    all_val = [current_val]

    for i in range(window_size, len(kmers_list)):
        current_val += kmers_value[kmers_list[i]] - kmers_value[kmers_list[i-window_size]]
        all_val.append(current_val)
    return np.array(all_val)/window_size

def window_slider_distance(kmers_list:list[int], kmers_ref_freq:Dict[int,float],  window_size:int=2000, p=2) -> np.array:
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


def KL_divergence(kmers_list:list[int], kmers_ref_freq:Dict[int,float], window_size:int=2000):
    """
    Computes Kullback-Leibler divergence between two probability distribution.
    Here computes the divergence between the window k-mer probability distribution and the overall 
    k-mer distribution in the genome.
    """
    init_kmer_window = kmers_list[0:window_size]
    window_freq = {kmer: init_kmer_window.count(kmer)/window_size for kmer in kmers_ref_freq.keys()}

    current_kl_div = np.sum(-kmers_ref_freq[kmer]*np.log10(freq/kmers_ref_freq[kmer]) if freq > 0 else 0 for kmer, freq in window_freq.items())
    all_kl_div = [current_kl_div]

    for i in range(window_size, len(kmers_list)):
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