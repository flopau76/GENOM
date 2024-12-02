from TP.loading import load_directory, load_directory_as_pointers
from TP.kmers import stream_kmers, stream_kmers_file
from time import time

def dict_intersection(dictA, dictB):
    """ Computes the intersection of two dictionaries
    :param dict dictA, dictB: dictionaries to compare"""
    intersection = 0
    for key, val in dictA:
        intersection += min(val, dictB.get(key, 0))
    return intersection

def list_intersection(listA, listB):
    """ Computes the intersection of two sorted lists
    :param np.array listA, listB: sorted np.array to compare"""
    intersection = 0
    idxA = 0
    idxB = 0
    while idxA < len(listA) and idxB < len(listB):
        if listA[idxA] == listB[idxB]:
            intersection += 1
            idxA += 1
            idxB += 1
        elif listA[idxA] < listB[idxB]:
            idxA += 1
        else:
            idxB += 1
    return intersection

def xorshift(val):
    """ Hash function using the xorshift algorithm """
    val ^= val << 13
    val &= 0xFFFFFFFFFFFFFFFF
    val ^= val >> 7
    val ^= val << 17
    val &= 0xFFFFFFFFFFFFFFFF
    return val


if __name__ == "__main1__":

    k = 21
    folder = "data"

    # Reading the data
    filenames = []
    kmers_list = []
    sequences = load_directory(folder)

    # Computing the kmers
    print("  Computing the kmers")
    for sample, seqs in sequences.items():
        print("Processing", sample)
        filenames.append(sample)
        kmers_list.append(sorted(list(stream_kmers(seqs, k))))

    # Computing the Jaccard index
    print("  Computing the pairwise similarities")
    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            intersection = list_intersection(kmers_list[i], kmers_list[j])
            dist_j = intersection / (len(kmers_list[i]) + len(kmers_list[j]) - intersection)
            print(filenames[i], filenames[j], dist_j)

if __name__ == "__main__":

    k = 21
    folder = "data"

    # Computing the kmers
    filenames = []
    kmers_list = []
    print("  Computing the kmers")
    for sample, file_pointer in load_directory_as_pointers(folder):
        print("Processing", sample)
        filenames.append(sample)
        kmers_list.append(sorted(list(stream_kmers_file(file_pointer, k))))

    # Computing the Jaccard index
    print("  Computing the pairwise similarities")
    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            intersection = list_intersection(kmers_list[i], kmers_list[j])
            dist_j = intersection / (len(kmers_list[i]) + len(kmers_list[j]) - intersection)
            print(filenames[i], filenames[j], dist_j)
