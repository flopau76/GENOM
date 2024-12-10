import numpy as np
from collections import Counter
from abc import ABC, abstractmethod

from typing import Dict, List, Tuple

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
    For it to work, we must have f(x1,...,xn) = h(g(x1), ..., g(xn)) 
    Having this, we can update f knowing only the changes in x1, ..., xn """

    @abstractmethod
    def compute_single(self, kmer:int, count:int) -> float:
        """ Compute the metric for a single kmer. """
        pass

    def compute(self, kmers_count:Dict[int,int]):
        """ Compute the metric for a window, assuming that if a kmer is not in kmers_count, it has no impact on the metric. """
        return sum([self.compute_single(kmer, count) for kmer, count in kmers_count.items()])

    def compute_initial(self, kmers_count:Dict[int,int]):
        """ Compute the metric for a window"""
        return self.compute(kmers_count)

    def update(self, old_value, count:Dict[int,int], changes:Dict[int,int]):
        """ Update the metric value after a window slide """
        old_count, new_count = get_subset(count, changes)
        return old_value + self.compute(new_count) - self.compute(old_count)

    def post_process(self, value:List[any]) -> np.ndarray:
        return np.array(value)
    
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
    
    def compute_initial(self, kmers_count:Dict[int,int]):
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
    
    def compute_initial(self, kmers_count:Counter[int]):
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
#                Iterate over file to compute metric                     #
##########################################################################

from io import TextIOWrapper
from TP.window_slider import stream_windows

def compute_metrics_file(file_pointer:TextIOWrapper, metric_list:List[Metric], window_size:int, k:int, step:int=1):
    """ Compute the metric for all windows of a file
    Currently supported metrics are: "average", "distance", "KLdiv", "jacc" """

    windows_stream = stream_windows(file_pointer, window_size, k, step=step)

    window, changes = next(windows_stream)
    results = [[metric.compute_initial(window.count) for metric in metric_list]]

    for window, changes in windows_stream:
        results.append([metric.update(old_value, window.count, changes) for metric, old_value in zip(metric_list, results[-1])])

    new_results = []
    for i, metric in enumerate(metric_list):
        res = [result[i] for result in results]
        new_results.append(metric.post_process(res))
    new_results = np.array(new_results)
    return new_results

##########################################################################
#       Old implementation using the kmer_list (for debugging only)      #
##########################################################################

def distance_kmers_list(kmers_list:List[int], ref_freq:Dict[int, float], window_size:int, p:int) -> np.array:
    kmers_count = Counter(kmers_list[:window_size])
    current_diff = {kmer: kmers_count[kmer]/window_size - ref for kmer, ref in ref_freq.items()}
    all_val = [sum([abs(diff)**p for diff in current_diff.values()])]
    for i in range(window_size, len(kmers_list)):
        current_diff[kmers_list[i]] = current_diff.get(kmers_list[i], 0) + 1/window_size
        current_diff[kmers_list[i-window_size]] -= 1/window_size
        all_val.append(sum([abs(diff)**p for diff in current_diff.values()]))
    return np.power(np.array(all_val),1/p)

def KLdivergence_kmers_list(kmers_list:List[int], kmers_freq:Dict[int, float], window_size:int) -> np.array:
    kmers_count = Counter(kmers_list[:window_size])
    current_freq = {kmer: count/window_size for kmer, count in kmers_count.items()}
    all_val = [sum([p*np.log10(p/q) for kmer, p in current_freq.items() if (q := kmers_freq.get(kmer, 0)) and p > 0])]
    for i in range(window_size, len(kmers_list)):
        current_freq[kmers_list[i]] = current_freq.get(kmers_list[i], 0) + 1/window_size
        current_freq[kmers_list[i-window_size]] -= 1/window_size
        all_val.append(sum([p*np.log10(p/q) for kmer, p in current_freq.items() if (q := kmers_freq.get(kmer, 0)) and p > 0]))
    return np.array(all_val)

def jaccard_kmers_list(kmers_list:List[int], ref_count:Dict[int, int], window_size:int) -> np.array:
    def jacc(countA, countB):
        return sum((countA & countB).values())/sum((countA | countB).values())
    current_count = Counter(kmers_list[:window_size])
    all_val = [jacc(current_count, ref_count)]
    for i in range(window_size, len(kmers_list)):
        current_count[kmers_list[i]] += 1
        current_count[kmers_list[i-window_size]] -= 1
        if current_count[kmers_list[i-window_size]] == 0:
            del current_count[kmers_list[i-window_size]]
        all_val.append(jacc(current_count, ref_count))
    return np.array(all_val)