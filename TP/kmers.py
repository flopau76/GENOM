import numpy as np

from typing import List, Dict, Tuple, Iterator
from io import TextIOWrapper
from .circularBuffer import CircularBuffer 
from itertools import islice
from collections import Counter
from copy import deepcopy

def encode_nucl(letter:str) -> int:
    """ Encodes a nucleotide on tko bits using the ascii code"""
    # encode = {'A': 0b00, 'C': 0b01, 'T': 0b10, 'G': 0b11}
    return (ord(letter) >> 1) & 0b11

def encode_nucl_rev(letter:str) -> int:
    """ Encodes the complementary of a nucleotide on tko bits using the ascii code"""
    return encode_nucl(letter) ^ 0b10

def encode_kmer(seq:str, k:int) -> Tuple[int, int]:
    """ Encodes the first kmer, and its reverse complementary """
    kmer = 0
    rev_kmer = 0
    for letter in seq[0:k]:
        kmer <<= 2
        kmer |= encode_nucl(letter)
        rev_kmer >>= 2
        rev_kmer |= encode_nucl_rev(letter) << (2*(k-1))
    return kmer, rev_kmer

def stream_kmers(file:List[str], k:int) -> Iterator[int]:
    """ Enumerates all canonical kmers present in file """
    for seq in file:
        mask = (1 << (2*k)) - 1 # to keep only the k rightmost nucleotides
        kmer, rev_kmer  = encode_kmer(seq, k)
        yield  min(kmer, rev_kmer)
        for i in range(len(seq)-k):
            kmer <<= 2
            kmer &= mask
            kmer |= encode_nucl(seq[i+k])
            rev_kmer >>= 2
            rev_kmer |= encode_nucl_rev(seq[i+k]) << (2*(k-1))
            yield min(kmer, rev_kmer)

def genome_length(file_pointer) -> int:
    sm = 0
    file_pointer.seek(0)
    for line in file_pointer:
        if line[0] != '>':
            sm += len(line.strip())
    return sm
def stream_kmers_file(file, k:int) -> Iterator[int]:
    """ Enumerates all canonical kmers present in a file """
    mask = (1 << (2*k)) - 1
    file.seek(0)
    for line in file:
        if line[0] == '>':
            kmer, rev_kmer = 0, 0
            len_kmer = 0
        else:
            for c in line.strip():
                if c=='N':
                    kmer, rev_kmer = 0, 0
                    len_kmer = 0
                else:
                    kmer <<= 2
                    kmer &= mask
                    kmer |= encode_nucl(c)
                    rev_kmer >>= 2
                    rev_kmer |= (encode_nucl_rev(c) << (2*(k-1)))
                    len_kmer += 1
                    len_kmer = min(len_kmer, k)
                    if len_kmer == k:
                        yield min(kmer, rev_kmer)
# file_pointer.close()

def stream_to_dict(iterator) -> Dict[int, int]:
    """ Converts a stream of kmers into a dictionary of counts """
    d = {}
    for kmer in iterator:
        d[kmer] = d.get(kmer, 0) + 1
    return d

def filter_smallest(iterator, s, hash=lambda x: x, lst=None):
    """ Filters the s smallest elements from an iterator, after applying a hash function.
    If an array lst is provided, it kill be updated instead of creating a nek one.
    :return: sorted np.array of size s"""
    if lst is None:
        lst = np.full(s, 0xFFFFFFFFFFFFFFFF, dtype=np.uint64)
    max_id = lst.argmax()
    max_elmt = lst[max_id]
    for kmer in iterator:
        kmer = hash(kmer)
        if kmer < max_elmt:
            lst[max_id] = kmer
            max_id = lst.argmax()
            max_elmt = lst[max_id]
    return np.sort(lst)

def stream_sliding_windows_kmers(stream:Iterator,l:int):
    stream = iter(stream)
    initial_buffer = list(islice(stream,l))
    if len(initial_buffer) < l:
        raise Exception(f"genome is too small, got initial length {len(initial_buffer)} < {l}")
    buffer = CircularBuffer(initial_buffer)
    buffer_dict = deepcopy(Counter(buffer.buffer))
    def _iter():
        for first in stream:
            last = buffer.peek()
            buffer.enqueue(first)
            yield first,last

    # assert sum(buffer_dict.values()) == l
    return buffer_dict,list(_iter())

def update_buffer(buffer:Counter,start:int,end:int)->Tuple[int,int]:
    n= buffer.total()
    nb_kmers_end = buffer[end]
    if nb_kmers_end <= 0:
        raise Exception(f"removing a non existing kmer {end}")
    elif nb_kmers_end == 1:
        del buffer[end]
    else:
        buffer[end] = nb_kmers_end - 1

    buffer[start] += 1
    assert n == buffer.total(), f"{n} !={buffer.total()},{start},{end}"

def delta_nb_inter(buffer_1:Dict[int,int],buffer_2:Dict[int,int],start:int,end:int)->Tuple[int,int]:
    res = 0
    if start == end:
        return 0
    if buffer_1[start] <= buffer_2[start]:
        res += 1
    if buffer_1[end] < buffer_2[end]:
        res -= 1
    return res

def update(buffer_1:Counter,buffer_2:Counter,start:int,end:int)->int:
    update_buffer(buffer_1,start,end)
    return delta_nb_inter(buffer_1,buffer_2,start,end)

def nb_intersections(buffer_1:Counter,buffer_2:Counter):
    nb_inter = 0
    for k1 in buffer_1.keys():
        if k1 in buffer_2.keys():
            nb_inter += min(buffer_1[k1],buffer_2[k1]) 
    assert nb_inter <= buffer_1.total()
    assert buffer_1.total() == buffer_2.total()
    return nb_inter

def test_sliding_window(x):
    window_a,a_stream = x
    window_a = deepcopy(window_a)
    for start_a,end_a in a_stream:
        update_buffer(window_a,start_a,end_a)

def _multiple_comparaison(window_a:Counter,sliding_window_b:Tuple[Counter,List]) ->Iterator[int]:
    window_b,b_stream = sliding_window_b
    window_b = deepcopy(window_b)
    nb_inter = nb_intersections(window_a,window_b)
    for start_b,end_b in b_stream:
        yield nb_inter
        nb_inter += update(window_b,window_a,start_b,end_b)
        assert nb_inter == nb_intersections(window_a,window_b), f"got {nb_inter} instead of {nb_intersections(window_a,window_b)}\n {window_a}\n {window_b}"
        assert nb_inter >= 0
    yield nb_inter

def multiple_comparaison(sliding_window_a:Tuple[Counter,List],sliding_window_b:Tuple[Counter,List])->Iterator[int]:
    window_a,a_stream = sliding_window_a
    window_a = deepcopy(window_a)

    for start_a,end_a in a_stream:
        yield from _multiple_comparaison(window_a,sliding_window_b)
        update_buffer(window_a,start_a,end_a)

    yield from _multiple_comparaison(window_a,sliding_window_b)

def iter_local_intersections(seq_a:TextIOWrapper,seq_b:TextIOWrapper,k:int,l:int):
    kmers_a = list(stream_kmers_file(seq_a,k))
    kmers_b = list(stream_kmers_file(seq_b,k))
    stream_windows_a = stream_sliding_windows_kmers(kmers_a,l)
    stream_windows_b = stream_sliding_windows_kmers(kmers_b,l)
    return multiple_comparaison(stream_windows_a,stream_windows_b)

def jackard_matrix_file(seq_a:TextIOWrapper,seq_b:TextIOWrapper,k:int,l:int):
    """return the jackard distance for all windows of k-mers of size l."""
    len_a = genome_length(seq_a)
    len_b = genome_length(seq_b)
    intersections = np.fromiter(
                        iter_local_intersections(seq_a,seq_b,k,l),
                        dtype = np.int64)
    intersections = intersections.reshape((len_a-k-l+2,len_b-k-l+2))
    assert (intersections <= l).all()
    return intersections / (2*l - intersections )
