from TP.loading import load_directory, load_directory_kmers
from TP.kmers import stream_kmers, filter_smallest
from time import time

def list_intersection(listA, listB):
    """ Computes the intersection of two sorted lists
    :param np.array listA, listB: sorted np.array to compare"""
    intersection = 0
    idxA = 0
    idxB = 0
    while idxA < listA.size and idxB < listB.size:
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


if __name__ == "__main__":

    k = 21
    s = 15000
    folder = "data"

    # Directly apply minhash when reading the data
    print("Reading the data and applying minhash")
    start = time()
    kmers = load_directory_kmers(folder, k, s, xorshift)
    print("  time:", time()-start)
    filenames = list(kmers.keys())
    kmers_lists = [kmers[filename] for filename in filenames]
    
    print("Computing the pairwise similarities")
    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            intersection = list_intersection(kmers_lists[i], kmers_lists[j])
            dist_j = intersection / (len(kmers_lists[i]) + len(kmers_lists[j]) - intersection)
            print(filenames[i], filenames[j], dist_j)
    print("  total time:", time()-start)