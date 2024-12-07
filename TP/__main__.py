from TP.loading import load_directory_as_pointers,iter_directory,open_genome
from TP.kmers import stream_kmers_file,jackard_matrix_file
import numpy as np
from typing import Dict, List
from time import time
from collections import Counter
import os, json
from itertools import product
import matplotlib.pyplot as plt


import TP.km_stats as km_stats
import TP.display as disp

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
        kmers_list = []
        dico = {}

        print("Processing", sample)
        size = len(file_pointer.read())
        file_pointer.seek(0)

        threshold = round(np.log2(size))
        
        kmers_list.append(list(stream_kmers_file(file_pointer, k)))

        kmer_split_list = np.array_split(*kmers_list, threshold)
        dico = dict(zip([f"{sample}_{i}" for i in range(threshold)], kmer_split_list))

        yield dico, kmers_list[0], sample


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
    folder = "toy_transfer"
    
    time_ = []
    out_dic = {}
    dir_path = os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]
    dir_path = '/'.join(dir_path)

    
    for dico_strain, kmers_list, sample in compute_kmer(folder, k):
        st = time()
        liste_jac = compute_jaccard(dico_strain)
        df = km_stats.load_as_matrix(liste_jac)

        print("  Starting frequence profile comparison")
        dico_km = Counter(kmers_list)
        freq_dico_km = {key:n/sum(list(dico_km.values())) for key, n in Counter(kmers_list).items()}

        hits, freq_avg_km = km_stats.window_slider(kmers_list, freq_dico_km)
        time_.append(time()-st)

        out_dic[sample] = hits
        #disp.display_freq(freq_avg_km)

    with open(os.path.join(dir_path,"transfer_summary.json"), 'w') as outjson:
        json.dump(out_dic, outjson)
    
    print("Average runtime per genome", np.average(time_[1:]))
    print("Total Runtime", round(sum(time()-st)))

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
    # for name_a,name_b in product(iter_directory(folder),iter_directory(folder)):
        # if name_a > name_b:
            # continue
        # print(f"comparing {name_a} and {name_b}")
    # with open_genome(name_a) as file_a,open_genome(name_b) as file_b:
        # m = jackard_matrix_file(file_a,file_b,k,l)
        # print(m)
