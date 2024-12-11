import numpy as np
from collections import Counter


class window:
    def __init__(self, kmers_stream):
        self._kmers = np.array(list(kmers_stream), dtype=int)
        self._head_position = 0
        self.size = len(self._kmers)
        self.count = Counter(self._kmers)

    def slide_right(self, new_kmer):
        old_kmer = self._kmers[self._head_position]
        self._kmers[self._head_position] = new_kmer
        self._head_position = (self._head_position + 1) % self.size
        self.count[new_kmer] += 1
        self.count[old_kmer] -= 1
        if self.count[old_kmer] == 0:
            del self.count[old_kmer]
        return old_kmer

    def __str__(self):
        return f"Window: {self._kmers[self._head_position]}...{self._kmers[(self._head_position-1)%self.size]}"