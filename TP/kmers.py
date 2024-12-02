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

def stream_slidinw_windows_kmers(stream:Iterator[int],l:int):
    buffer = circular_buffer(first(l,stream))
    yield circular_buffer_as_dict(buffer)
    for first in buffer:
        last = buffer.circle(first)
        yield first,last

def update_start(buffer_1:dict,buffer_2:dict,start:int)->Tuple[int,int]:
    n1 =  buffer_1.get(start,0)
    buffer_1[start] = n1 + 1
    if n1 >= buffer_2[start]:
        return 1,0
    else:
        return 0,1

def update_end(buffer_1:Dict[int,int],buffer_2:Dict[int,int],end:int)->Tuple[int,int]:
    n1 =  buffer_1.get(end,0)
    if n1 == 0:
        error(f"removing a non existing kmer {end}")
    elif n1 == 1:
        del buffer[end]
    else:
        buffer_1[end] = n1 - 1

    if n1 > buffer_2[start]:
        return -1,0
    else:
        return 0,-1

def update(buffer_1:Dict[int,int],buffer_2:Dict[int,int],start:int,end:int)->Tuple[int,int]:
    delta_union_1,delta_inter_1  = update_start(buffer_1,buffer_2,start)
    delta_union_2,delta_inter_2  = update_end(buffer_1,buffer_2,start)
    return delta_union_1 + delta_union_2 ,delta_inter_1 + delta_inter_2

def nb_union_inter(buffer_1:Dict[int,int],buffer_2:Dict[int,int],):
    nb_union,nb_inter = 0,0
    for k1 in buffer_1.keys():
        if k1 in buffer_2.keys():
            nb_union += max(buffer_1[k1],buffer_2[k1]) 
            nb_inter += min(buffer_1[k1],buffer_2[k1]) 
        else:
            nb_union += buffer_1[k1]
            nb_inter += buffer_1[k1]
    for k2 in buffer_2.keys():
        if k2 not in buffer_1.keys():
            nb_union += buffer_1[k1]
            nb_inter += buffer_1[k1]
    return nb_union,nb_inter

def multiple_comparaison(a_stream:Iterator,b_stream:Iterator):
    nb_union,nb_inter = nb_union_inter(a_start,b_start)
    a_buffer = next(a_stream)
    b_buffer = next(b_stream)
    yield (nb_union,nb_inter)
    for a_s,a_e in a_stream:
        delta_union, delta_inter = update(a_buffer,b_buffer,a_s,a_e)
        nb_union += delta_union
        nb_inter += delta_inter
        for b_s,b_e in b_stream:
            delta_union,delta_inter = update(b_buffer,a_buffer,b_s,b_e)
            nb_union += delta_union
            nb_inter += delta_inter
            yield (nb_union,nb_inter)
