import gzip
from os import listdir, path

from TP.kmers import filter_smallest, stream_kmers_file


def load_fasta(file_pointer):
    """ Loads a fasta formated file into a list of sequences.
    :param file_pointer: The stream of the fasta file to load.
    :return Array: An array of strings where each string is a sequence from the fasta
    """
    texts = []
    txt = []

    for line in file_pointer:
        if line[0] == '>':
            if len(txt) > 0:
                texts.append("".join(txt))
            txt = []
        else:
            txt.append(line.strip())

    if len(txt) > 0:
        texts.append("".join(txt))
    return texts


def load_directory(directory):
    """ Loads all the fasta files from a data directory into a dictionary.
    Each subdirectory in data is considered as a different sample.
    Fatsta files (even gzipped) are loaded.
    :param str directory: Path to the data directory to load.
    :return dict: A dict containing pairs (sample, [sequence list]).
    """
    sequence_dict = {}
    for name in listdir(directory):
        subpath = path.join(directory, name)
        # Look for sample directories
        if path.isdir(subpath):
            # Creates one list of sequence per sample
            sequence_dict[name] = []
            for filename in listdir(subpath):
                # Load raw fasta files
                if filename.endswith(".fa") or filename.endswith(".fasta"):
                    with open(path.join(subpath, filename)) as fp:
                        sequence_dict[name] += load_fasta(fp)
                        print("Loaded", filename, len(sequence_dict[name]))
                # Load gzipped fasta files
                elif filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
                    with gzip.open(path.join(subpath, filename), 'rt') as fp:
                        sequence_dict[name] += load_fasta(fp)
                        print("Loaded", filename, len(sequence_dict[name]))
    
    return sequence_dict


def load_directory_kmers(directory, k, s, hash=lambda x: x):
    """ Loads all the fasta files from a data directory into a dictionary.
    Instead of keeping all sequences in memory, we directly compute the kmers and minhash them
    :param str directory: Path to the data directory to load.
    :param int k: Size of kmers to compute
    :param int s: Number of kmers to keep per sample
    :param function hash: Hash function to use on kmers
    :return dict: A dict containing pairs (sample, [minhashed kmer list]).
    """
    kmers_dict = {}
    for name in listdir(directory):
        subpath = path.join(directory, name)
        # Look for sample directories
        if path.isdir(subpath):
            # Creates one list of sequence per sample
            kmers_lst = None
            for filename in listdir(subpath):
                # Load raw fasta files
                file_pointer = None
                if filename.endswith(".fa") or filename.endswith(".fasta"):
                    file_pointer = open(path.join(subpath, filename))
                # Load gzipped fasta files
                elif filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
                    file_pointer = gzip.open(path.join(subpath, filename), 'rt')
                if file_pointer is not None:
                    kmers_lst = filter_smallest(stream_kmers_file(file_pointer, k), s, hash, kmers_lst)
                    file_pointer.close()
            kmers_dict[name] = kmers_lst
            print("Loaded", name)
    
    return kmers_dict


if __name__ == "__main__":
    files = load_directory("data")
    print(len(files))
