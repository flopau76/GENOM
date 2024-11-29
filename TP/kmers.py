import numpy as np

from typing import List, Dict, Tuple, Iterator
from io import TextIOWrapper

def encode_nucl(letter:str) -> int:
    """ Encodes a nucleotide on two bits using the ascii code"""
    # encode = {'A': 0b00, 'C': 0b01, 'T': 0b10, 'G': 0b11}
    return (ord(letter) >> 1) & 0b11

def encode_nucl_rev(letter:str) -> int:
    """ Encodes the complementary of a nucleotide on two bits using the ascii code"""
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

def stream_kmers_file(file_pointer:TextIOWrapper, k:int) -> Iterator[int]:
    """ Enumerates all canonical kmers present in a file """
    mask = (1 << (2*k)) - 1
    for line in file_pointer:
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
    file_pointer.close()

def stream_to_dict(iterator) -> Dict[int, int]:
    """ Converts a stream of kmers into a dictionary of counts """
    d = {}
    for kmer in iterator:
        d[kmer] = d.get(kmer, 0) + 1
    return d

def filter_smallest(iterator, s, hash=lambda x: x, lst=None):
    """ Filters the s smallest elements from an iterator, after applying a hash function.
    If an array lst is provided, it will be updated instead of creating a new one.
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

def multiple_comparaison(a_stream:Iterator,b_stream:Iterator,a_start,b_start):
    nb_union,nb_inter = nb_union_inter(a_start,b_start)
    a_buffer = copy(a_start)
    b_buffer = copy(b_start)
    yield (nb_union,nb_inter)
    for a_s,a_e in a_stream:
        nb_union += update_union_s(a_buffer,b_buffer,a_s)
        nb_inter += update_inter_s(a_buffer,b_buffer,a_s)
        nb_union += update_union_e(a_buffer,b_buffer,a_s)
        nb_inter += update_inter_e(a_buffer,b_buffer,a_s)
        for b_s,b_e in b_stream:
            nb_union += update_union_s(b_buffer,a_buffer,a_s)
            nb_inter += update_inter_s(b_buffer,a_buffer,a_s)
            nb_union += update_union_e(b_buffer,a_buffer,a_s)
            nb_inter += update_inter_e(b_buffer,a_buffer,a_s)
            yield (nb_union,nb_inter)
