from TP.loading import load_directory_as_pointers,iter_directory,open_genome
from TP.kmers import stream_kmers_file,jackard_matrix_file
import numpy as np
from typing import Dict, List
from time import time
from collections import Counter
import os, json
from itertools import product
import matplotlib.pyplot as plt
import TP.signatures as signatures
import TP.km_stats as km_stats
import TP.display as disp

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.signal import find_peaks
from time import time
import os, json

from typing import Dict, List

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

def compute_kmer(folder : str, k : int):
    # Computing the kmers
    print("  Computing the kmers")
    for sample, file_pointer in load_directory_as_pointers(folder):
        print("Processing", sample)        
        yield sample, list(stream_kmers_file(file_pointer, k))

def compute_jaccard(dico : Dict[str, List[int]]):
    # Computing the Jaccard index
    print("  Computing the pairwise similarities")
    filenames = list(dico.keys())
    list_tuple_jac = []

    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):
            intersection = list_intersection(sorted(dico[filenames[i]].tolist()), sorted(dico[filenames[j]].tolist()),k=8,l=16)
            dist_j = intersection / (len(dico[filenames[i]]) + len(dico[filenames[j]]) - intersection)
            #print(f"{'==='*20}\n{filenames[i]} | {filenames[j]} | {dist_j}")

            list_tuple_jac.append((filenames[i], filenames[j], dist_j))

    return list_tuple_jac


def dump_matrix_to_csv(matrix, filename, delimiter=',', precision=6):
    """
    Dump a NumPy matrix to a CSV file.
    
    Parameters:
    -----------
    matrix : numpy.ndarray
        The input matrix to dump
    filename : str 
        Path to save the CSV file. 
    delimiter : str, optional
        Delimiter to use in the CSV file (default: ',')
    precision : int, optional
        Number of decimal places to use when writing float values (default: 6)
    
    Returns:
    --------
    str
        Path to the saved CSV file
    """
    
    # Ensure the filename has .csv extension
    if not filename.lower().endswith('.csv'):
        filename += '.csv'
    
    # Use numpy's savetxt for straightforward CSV export
    try:
        # Set print options to control precision
        np.set_printoptions(precision=precision, suppress=True)
        
        # Save the matrix to CSV
        np.savetxt(filename, matrix, delimiter=delimiter, fmt=f'%0.{precision}f')
        
        print(f"Matrix successfully exported to: {filename}")
        
        # Additional matrix information
    
    except Exception as e:
        print(f"Error exporting matrix to CSV: {e}")

def visualize_matrix(matrix, 
                    save_path, 
                    title='Matrix Heatmap', 
                    cmap='viridis', 
                    figsize=(10, 10), 
                    value_fmt='%.2f',
                    colorbar=True,
                    text_color='black'):
    """
    Visualize a NumPy matrix as a heatmap using Matplotlib.
    
    Parameters:
    -----------
    matrix : numpy.ndarray
        Input matrix to visualize
    title : str, optional
        Title of the heatmap (default: 'Matrix Heatmap')
    cmap : str, optional
        Colormap to use (default: 'viridis')
    save_path : str 
        Path to save the heatmap image (default: None)
    figsize : tuple, optional
        Figure size (width, height) in inches (default: (10, 10))
    value_fmt : str, optional
        Format specifier for cell values (default: '%.2f')
    colorbar : bool, optional
        Whether to show a colorbar (default: True)
    text_color : str, optional
        Color of the text in cells (default: 'black')
    
    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure object
    """
    # Create a new figure with specified size
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the heatmap
    im = ax.imshow(matrix, cmap=cmap, aspect='equal',interpolation ='none')
    
    # Set title
    ax.set_title(title)
    
    # Add colorbar
    if colorbar:
        plt.colorbar(im, ax=ax, label='Value')
    
    
    # Save figure if path is provided
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    
    return fig

if __name__ == "__main__1":
    k = 8
    window_size = 2000
    input_folder = "toy_no_transfer"
    output_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", f"transfer_summary_{input_folder}.json")

    best_hits = {}

    for sample, file_pointer in  load_directory_as_pointers(input_folder):
        start = time()
        print("Processing sample ", sample)
        print("    Computing the list of kmers")
        kmers_list = list(stream_kmers_file(file_pointer, k))

        print("    Starting frequence profile comparison")
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}
        
        st = time()
        window_distance = signatures.window_slider_distance(kmers_list, kmers_freq, window_size=window_size)
        print("   L2 done in ", round(time()-st, 4), "s")

        st = time()
        window_kldiv = signatures.KLdivergence(kmers_list, kmers_freq)
        print("   KL divergencce done in ", round(time()-st, 4), "s")
        disp.display_windows(window_kldiv, ylabel="KL divergence", title="KL divergence for sample "+sample)
        
        #kldiv2 = signatures.naive_KLdiv(kmers_list, kmers_freq)
        #disp.display_freq(kldiv2)
        
        # kmers_rarity = {kmer:1/freq for kmer, freq in kmers_freq.items()}
        # window_average_rarity = signatures.window_slider_average(kmers_list, kmers_rarity, window_size)

        print("    Finding the best hits")
        # highest_values_indices = signatures.find_maxima(window_distance, nb_hits)
        # highest_values = window_distance[highest_values_indices]

        """
        # possible filtering parameters: dependig on the metric, they may need to be adjusted
        height = np.mean(window_distance) + 10*np.var(window_distance)
        distance = window_size"""
        prominence = (np.max(window_distance) - np.min(window_distance))/3
        highest_values_indices, _ = find_peaks(window_distance, prominence=prominence)
        highest_values = window_distance[highest_values_indices]
        best_hits[sample] = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

        print("    Done in ", round(time()-start, 4), "s")

    plt.waitforbuttonpress()
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)


if __name__ == "__main__":
    k = 8
    l = 16
    folder = "toy_transfer"
    out= "result"
    if not(os.path.exists(out)):
        os.mkdir(out)
    
    dir_path = os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]
    dir_path = '/'.join(dir_path)

    
    for name_a in iter_directory(folder):
        with open_genome(name_a) as file_a:
            out_path= os.path.join(out,os.path.splitext(os.path.basename(name_a))[0])
            m = jackard_matrix_file(file_a,file_a,k,l)
            dump_matrix_to_csv(m,out_path)
            visualize_matrix(m,out_path)
            print(m)
    for name_a,name_b in product(iter_directory(folder),iter_directory(folder)):
        if name_a >= name_b:
            continue
        print(f"comparing {name_a} and {name_b}")
        with open_genome(name_a) as file_a,open_genome(name_b) as file_b:
            out_path= os.path.join(out,
                os.path.splitext(os.path.basename(name_a))[0] +
                "_" +
                os.path.splitext(os.path.basename(name_b))[0])
            m = jackard_matrix_file(file_a,file_b,k,l)
            dump_matrix_to_csv(m,out_path)
            visualize_matrix(m,out_path)
            print(m)
