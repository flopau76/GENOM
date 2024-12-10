import numpy as np
from collections import Counter
from itertools import islice

from TP.kmers import stream_kmers_file

from typing import Dict, List, Iterator, Tuple, Set
from io import TextIOWrapper

class window:
    def __init__(self, kmers_stream):
        self.position = 0
        self.kmers = np.array(list(kmers_stream))
        self._head_position = 0
        self.count = Counter(self.kmers)
        self.size = len(self.kmers)

    def slide_right(self, new_kmer):
        self.position += 1
        old_kmer = self.kmers[self._head_position]
        self.kmers[self._head_position] = new_kmer
        self._head_position = (self._head_position + 1) % self.size
        self.count[new_kmer] += 1
        self.count[old_kmer] -= 1
        if self.count[old_kmer] == 0:
            del self.count[old_kmer]
        return old_kmer

    def __str__(self):
        return f"Window at position {self.position}: {self.kmers[self._head_position]}...{self.kmers[(self._head_position-1)%self.size]}"

def stream_windows(file_pointer:TextIOWrapper, window_size:int, k:int, step:int=1) -> Iterator[Tuple[window, Counter, Counter]]:
    """ Enumerate all windows of window_size kmers in a file, and the kmers that changed since the last window
    If step is larger than 1, only every step windows are returned """
    kmers_stream = stream_kmers_file(file_pointer, k)
    current_window = window(islice(kmers_stream, window_size))
    assert len(current_window.kmers) == window_size, f"Initial window size {window_size} is larger than genome size {len(current_window.kmers)}"
    yield current_window, None

    count = 0
    changes = dict()
    for new_kmer in kmers_stream:
        old_kmer = current_window.slide_right(new_kmer)
        changes[old_kmer] = changes.get(old_kmer, 0) - 1
        changes[new_kmer] = changes.get(new_kmer, 0) + 1
        count += 1
        if count % step == 0:
            yield current_window, changes
            count = 0
            changes = dict()