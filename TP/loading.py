import gzip
from os import listdir, path

from typing import Iterator, Tuple, List, Dict
from io import TextIOWrapper


def load_fasta(file_pointer:TextIOWrapper) -> List[str]:
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


def load_directory(directory:str) -> Dict[str, List[str]]:
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

def open_genome(name:str):
    if name.endswith(".fa") or name.endswith(".fasta"):
        return open(name)
    elif name.endswith(".fa.gz") or name.endswith(".fasta.gz"):
        return gzip.open(name, 'rt')
    else:
        raise Exception(f"file {name} has wrong format")

def iter_directory(directory:str):
    for name in listdir(directory):
        subpath = path.join(directory, name)
        if path.isdir(subpath):
            for filename in listdir(subpath):
                yield path.join(subpath,filename)
    
def load_directory_as_pointers(directory:str) -> Iterator[Tuple[str, TextIOWrapper]]:
    """ Creates a generator that yields the name of the sample and a file pointer to the fasta file. """
    for name in listdir(directory):
        subpath = path.join(directory, name)
        if path.isdir(subpath):
            for filename in listdir(subpath):
                if filename.endswith(".fa") or filename.endswith(".fasta"):
                    yield (name, open(path.join(subpath, filename)))
                elif filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
                    yield (name, gzip.open(path.join(subpath, filename), 'rt'))


if __name__ == "__main__":
    files = load_directory("data")
    print(len(files))
