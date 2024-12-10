import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, deque
from itertools import islice
from time import time
import os, json

from scipy.signal import find_peaks
from scipy.special import kl_div

from typing import Dict, List, Iterator, Tuple
from io import TextIOWrapper

from TP.loading import iter_directory, open_genome
from TP.kmers import stream_kmers_file
from TP.display import display_windows

class window:
    def __init__(self, kmers_stream):
        self.kmers = deque()
        self.count = Counter()
        self.position = 0
        self.size = 0
        for kmer in kmers_stream:
            self.kmers.append(kmer)
            self.count[kmer] += 1
            self.size += 1
        self.old_kmer = None
        self.new_kmer = None

    def slide_right(self, new_kmer):
        self.new_kmer = new_kmer
        self.kmers.append(self.new_kmer)
        self.old_kmer = self.kmers.popleft()
        self.count[self.new_kmer] += 1
        self.count[self.old_kmer] -= 1
        if self.count[self.old_kmer] == 0:
            del self.count[self.old_kmer]
        self.position += 1

    def __str__(self):
        return f"Window at position {self.position}, size {self.size}, old kmer {self.old_kmer}, new kmer {self.new_kmer}"

def stream_windows(file_pointer:TextIOWrapper, window_size:int, k:int, step:int=1) -> Iterator[window]:
    """ Enumerate all windows of window_size kmers in a file
    If step is larger than 1, only every step windows are returned """
    kmers_stream = stream_kmers_file(file_pointer, k)
    current_window = window(islice(kmers_stream, window_size))
    assert len(current_window.kmers) == window_size, f"Initial window size {window_size} is larger than genome size {len(current_window.kmers)}"
    yield current_window

    count = 0
    for new_kmer in kmers_stream:
        current_window.slide_right(new_kmer)
        count += 1
        if count % step == 0:
            count = 0
            yield current_window

def compute_average(window:window, ref_value:Dict[int, float]) -> float:
    """ Compute the average value of the kmers in the window """
    return sum([ref_value[kmer] for kmer in window.kmers])/window.size

def update_average(old_average:float, new_window:window, ref_value:Dict[int, float]) -> np.array:
    """ Efficiently update the average value of the window """
    return old_average + (ref_value[new_window.new_kmer] - ref_value[new_window.old_kmer])/new_window.size

def compute_distance(window:window, ref_freq:Dict[int, float], p=2) -> float:
    """ Compute the average Lp distance of the kmers frequencies in the window, to given reference kmer frequencies 
    Note: for efficient convolution, we return d**p instead of d """
    diff = [ window.count[kmer]/window.size - ref_freq.get(kmer,0) for kmer in set(window.count.keys()).union(ref_freq.keys())]
    diff = np.array(diff)
    return sum(np.power(np.abs(diff),p))

def update_distance(old_distance:float, new_window:window, ref_freq:Dict[int, float], p=2) -> float:
    """ Efficiently update the Lp distance of the window """
    def distance(a,b):
        return np.power(np.abs(a-b),p)
    if new_window.old_kmer == new_window.new_kmer:
        new_distance = old_distance
    else:
        new_distance = old_distance \
                  + distance(new_window.count[new_window.new_kmer]/new_window.size, ref_freq[new_window.new_kmer]) \
                  + distance(new_window.count[new_window.old_kmer]/new_window.size, ref_freq[new_window.old_kmer]) \
                  - distance((new_window.count[new_window.new_kmer]-1)/new_window.size, ref_freq[new_window.new_kmer]) \
                  - distance((new_window.count[new_window.old_kmer]+1)/new_window.size, ref_freq[new_window.old_kmer])
    return new_distance

def post_process_distance(distance:float, p:int=2) -> float:
    return np.power(distance, 1/p)

def compute_KLdivergence(window:window, ref_freq:Dict[int, float]) -> float:
    """ Compute the Kullback-Leibler divergence of the kmers frequencies in the window, to given reference kmer frequencies """
    count = 0
    for kmer, p in window.count.items():
        p = p/window.size
        q = ref_freq[kmer]    # note: if the window has a kmer not in the ref_freq (q=0), the KLdivergence will be infinite
        count += p*np.log10(p/q)
    return count

def update_KLdivergence(old_KLdiv:float, new_window:window, ref_freq:Dict[int, float]) -> float:
    """ Efficiently update the Kullback-Leibler divergence of the window """
    def KLdiv(p,q):
        return p*np.log10(p/q) if p > 0 else 0
    if new_window.old_kmer == new_window.new_kmer:
        new_KLdiw = old_KLdiv
    else:
        new_KLdiw = old_KLdiv \
                + KLdiv(new_window.count[new_window.new_kmer]/new_window.size, ref_freq[new_window.new_kmer]) \
                + KLdiv(new_window.count[new_window.old_kmer]/new_window.size, ref_freq[new_window.old_kmer]) \
                - KLdiv((new_window.count[new_window.new_kmer]-1)/new_window.size, ref_freq[new_window.new_kmer]) \
                - KLdiv((new_window.count[new_window.old_kmer]+1)/new_window.size, ref_freq[new_window.old_kmer])
    return new_KLdiw

def compute_jaccard(window:window, ref_count:Dict[int, int]) -> float:
    """ Compute the Jaccard index of the kmers in the window, to given reference kmer counts 
    Note: for efficient convolution, we return the cardinal of the intersection and union instead of their ratio"""
    inter = 0
    union = 0
    kmer_set = set(window.count.keys()).union(ref_count.keys())
    for kmer in kmer_set:
        inter += min(window.count[kmer], ref_count.get(kmer,0))
        union += max(window.count[kmer], ref_count.get(kmer,0))
    return (inter, union)

def update_jaccard(old_jaccard:Tuple[int,int], new_window:window, ref_count:Dict[int, int]) -> float:
    """ Efficiently update the Jaccard index of the window """
    if new_window.old_kmer == new_window.new_kmer:
        new_jaccard = old_jaccard
    else:
        inter, union = old_jaccard
        if new_window.count[new_window.new_kmer] > ref_count[new_window.new_kmer]:
            union += 1
        else:
            inter += 1
        if new_window.count[new_window.old_kmer] < ref_count[new_window.old_kmer]:
            inter -= 1
        else:
            union -= 1
        new_jaccard = (inter, union)
    return new_jaccard

def post_process_jaccard(jaccard:Tuple[int,int]) -> float:
    """ Compute the Jaccard index from the cardinal of the intersection and union """
    jaccard = np.array(jaccard)
    return jaccard[:,0]/jaccard[:,1]

def compute_metrics_file(file_pointer:TextIOWrapper, metrics:List[str], window_size:int, k:int, step:int=1, **kwargs):
    """ Compute the metric for all windows of a file
    Currently supported metrics are: "average", "distance", "KLdiv", "jacc" """
    def compute_metric(window:window):
        value = []
        for met in metrics:
            if met == "average":
                assert "ref_value" in kwargs, "Reference value is needed for average metric"
                value.append(compute_average(window, kwargs["ref_value"]))
            elif met == "distance":
                assert "ref_freq" in kwargs, "Reference frequency is needed for distance metric"
                p = kwargs.get("p", 2)
                value.append(compute_distance(window, kwargs["ref_freq"], p=p))
            elif met == "KLdiv":
                assert "ref_freq" in kwargs, "Reference frequency is needed for KLdivergence metric"
                value.append(compute_KLdivergence(window, kwargs["ref_freq"]))
            elif met == "jacc":
                assert "ref_count" in kwargs, "Reference count is needed for Jaccard metric"
                value.append(compute_jaccard(window, kwargs["ref_count"]))
            else:
                raise ValueError(f"Unknown metric {met}")
        return value

    def update_metric(old_val:List[float], new_window:window):
        new_val = []
        for i, met in enumerate(metrics):
            if met == "average":
                new_val.append(update_average(old_val[i], new_window, kwargs["ref_value"]))
            elif met == "distance":
                p = kwargs.get("p", 2)
                new_val.append(update_distance(old_val[i], new_window, kwargs["ref_freq"], p=p))
            elif met == "KLdiv":
                new_val.append(update_KLdivergence(old_val[i], new_window, kwargs["ref_freq"]))
            elif met == "jacc":
                new_val.append(update_jaccard(old_val[i], new_window, kwargs["ref_count"]))
        return new_val
    
    def post_process_metric(result:List[any]):
        new_result = []
        for i, met in enumerate(metrics):
            res = [res[i] for res in result]
            if met == "distance":
                p = kwargs.get("p", 2)
                new_result.append(post_process_distance(res, p))
            elif met == "jacc":
                new_result.append(post_process_jaccard(res))
            else:
                new_result.append(res)
        return np.array(new_result).T

    windows_stream = stream_windows(file_pointer, window_size, k, step=step)

    new_window = next(windows_stream)
    results = [compute_metric(new_window)]
    for new_window in windows_stream:
        results.append(update_metric(results[-1], new_window))

    return post_process_metric(results)


def find_maxima(window_array:np.array, nb_maxima:int) -> np.array:
    """Finds the indices of the nb_maxima highest values in the window_array"""
    return np.argpartition(-window_array, nb_maxima)[:nb_maxima]

def find_local_maxima(window_array:np.array, **kwargs) -> np.array:
    """Finds the indices of the local maxima above a given threshold in the window_array"""
    return find_peaks(window_array, **kwargs)[0]


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

##########################################################################

if __name__ == "__main__":
    k = 8
    window_size = 2000
    input_folder = "toy_transfer"
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    best_hits = {}

    metrics = ["distance", "KLdiv", "jacc"]

    for sample, file_name in  iter_directory(input_folder):
        file_pointer = open_genome(file_name)

        # iterate once over the file to compute the total kmer frequency
        kmers_list = list(stream_kmers_file(file_pointer, k))
        kmers_count = Counter(kmers_list)
        ref_count = Counter(kmers_list[:window_size])
        kmers_nb = sum(kmers_count.values())
        kmers_freq = {kmer:count/kmers_nb for kmer, count in kmers_count.items()}

        results = compute_metrics_file(file_pointer, metrics, window_size, k, ref_freq=kmers_freq, ref_count=ref_count)
        file_pointer.close()

        fig = plt.figure()
        for i, metric in enumerate(metrics):
            ax = fig.add_subplot(len(metrics), 1, i+1)
            display_windows(results[:,i], title=metrics[i], ax=ax)
            ax.set_xlabel(None)
        ax.set_xlabel("Position")
        
        fig.suptitle(f"Metrics for sample {sample}")
        fig.tight_layout()

        fig.savefig(f"transfer_summary_{sample}_test.png")

        # storing best hits:
        results = results[:,0] # distance
        highest_values_indices, _ = find_peaks(results, prominence=np.max(results) - np.min(results)/3)
        highest_values = results[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

    # dumping hits to json
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)

    # wainting for all figures to be closed manually
    plt.show(block=True)