from typing import List, Dict, Tuple, Iterator
from io import TextIOWrapper

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

def genome_length(file_pointer:TextIOWrapper) -> int:
    """ Returns the length of the genome in a fasta file """
    sm = 0
    file_pointer.seek(0)
    for line in file_pointer:
        if line[0] != '>':
            sm += len(line.strip())
    return sm

def stream_kmers_file(file_pointer:TextIOWrapper, k:int) -> Iterator[int]:
    """ Enumerates all canonical kmers present in a file """
    file_pointer.seek(0)
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