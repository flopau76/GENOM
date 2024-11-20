import numpy as np

def encode_nucl(letter):
    """ Encodes a nucleotide on two bits
    using the ascii code"""
    # encode = {'A': 0b00, 'C': 0b01, 'T': 0b10, 'G': 0b11}
    return (ord(letter) >> 1) & 0b11

def encode_nucl_rev(letter):
    """ Encodes the complementary of a nucleotide on two bits
    using the ascii code"""
    return encode_nucl(letter) ^ 0b10

def encode_kmer(seq, k):
    """ Encodes the first kmer, and its reverse complementary """
    kmer = 0
    rev_kmer = 0
    for letter in seq[0:k]:
        kmer <<= 2
        kmer |= encode_nucl(letter)
        rev_kmer >>= 2
        rev_kmer |= encode_nucl_rev(letter) << (2*(k-1))
    return kmer, rev_kmer

def stream_kmers(file, k):
    """ Enumerates all kmers present in file (list of sequences)
    Note: we keep only the lowest value of each kmer and its reverse complementary """
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

def stream_kmers_file(file_pointer, k):
    """ Enumerates all canonical kmers present in a fasta file
    :param file_pointer: The stream of the fasta file to load.
    """
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