import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
from scipy.signal import find_peaks

from collections import Counter
import argparse
import os, json, time
from typing import Dict

import compute_signatures.metrics as metrics
from compute_signatures.loading import iter_directory, open_genome
from compute_signatures.kmers import stream_kmers_file
import compute_signatures.display as display
from compute_signatures.ribo_sort import drop_ribo_position


def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'
    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    if iteration == total: 
        print()

def find_potential_HGT(result:np.ndarray, min_height=None, min_prominence=None, path_ribo_db=None, sample=None, ribo_genome_file_table:dict={}) -> Dict[int, float]:
    """ Find potential HGT regions corresponding to the highest peaks in the result """
    highest_values_indices, _ = find_peaks(result, prominence=min_prominence, height=min_height)
    highest_values = result[highest_values_indices]
    best_hits = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}

    if path_ribo_db is not None:
        ribo_file_name = ribo_genome_file_table[sample]
        ribo_file_path = os.path.join(path_ribo_db, ribo_file_name)
        best_hits = drop_ribo_position(best_hits, ribo_file_path)

    return best_hits


if __name__ == "__main__":
    metric_dict = {0:metrics.distance(), 1:metrics.chi_squared(), 2:metrics.KLdivergence(), 3:metrics.Convolution()}

    parser = argparse.ArgumentParser(description="Compute the signature of a genome and find potential HGT regions")
    parser.add_argument("input_db", help="The name of the input database (must be in `input/`)")
    parser.add_argument('-o', '--output', help="The name of the output (will be in `output/transfer_summary`)", type=str, default=None)
    parser.add_argument('-k', '--kmer', help='The size of the kmer (default=8)', type=int, default=5)
    parser.add_argument('-w', '--window', help='The size of the sliding window (default=2000)', type=int, default=5000)
    parser.add_argument('-r', '--ref', help='Path to the file containing the known HGT for them to be shown on the report', type=str, default=None)
    parser.add_argument('-m', '--metric', help=f'Metric used for computation. Currently supports: ' + " ,".join([f"{metric.name} ({key})" for key, metric in metric_dict.items()]), type=int, default=1)
    parser.add_argument('-b', '--ribo', help="Path to the ribosome database if you wish to filter them out", type=str, default=None)
    
    args = parser.parse_args()

    k = args.kmer
    window_size = args.window
    input_name = args.input_db
    reference_dico = args.ref
    if reference_dico is not None:
        reference_dico = display.ref_parse(reference_dico)

    output_name = args.output
    metric = metric_dict[args.metric]
    path_ribo_db = args.ribo

    if output_name is None:
        output_name = os.path.basename(input_name) + '_' + metric.name

    base_dir =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    input_folder = os.path.join(base_dir, "input", input_name)
    output_path_json = os.path.join(base_dir, "output", f"transfer_summary", f"{output_name}.json")
    output_path_pdf = os.path.join(base_dir, "output", f"transfer_summary", f"{output_name}.pdf")

    ribo_genome_file_table = None
    if path_ribo_db is not None:
        liste_dir = os.listdir(input_folder)
        ribo_genome_file_table = dict(zip(liste_dir, os.listdir(path_ribo_db)))

    pdf = bpdf.PdfPages(output_path_pdf)

    best_hits = {}

    times_kmers = []
    times_windows = []
    start = time.time()

    files = iter_directory(input_folder)
    n = 0
    n_total = len(files)
    progressbar(n, n_total , prefix = 'Progress:', suffix = 'Complete', printEnd = "\r")

    for sample, file_name in files:
        # iterate once over the file to get the kmers_list
        t0 = time.time()
        with open_genome(file_name) as file_pointer:
            kmers_list = list(stream_kmers_file(file_pointer, k))
        
        # compute the average signature of the genome
        kmers_count = Counter(kmers_list)
        kmers_freq = {kmer:count/len(kmers_list) for kmer, count in kmers_count.items()}
        t1 = time.time()
        times_kmers.append(t1-t0)

        # compute the distance to the average signature along the genome
        metric.update(ref_freq=kmers_freq, ref_count=kmers_count)
        result = metric.slide_window(kmers_list, window_size)
        t2 = time.time()
        times_windows.append(t2-t1)

        # find the highest peaks
        sample_hits = find_potential_HGT(result, min_prominence=(np.max(result) - np.min(result))/3, path_ribo_db=path_ribo_db, sample=sample, ribo_genome_file_table=ribo_genome_file_table)
        best_hits[sample] = sample_hits

        # save the resulting figure
        fig = display.display_windows(result, sample, hits=sample_hits, title=f"{sample}", ylabel=metric.name, dpi=300, ref=reference_dico, window_size=window_size)
        fig.savefig(pdf, format='pdf')
        plt.close(fig)

        n += 1
        progressbar(n, n_total, prefix = 'Progress:', suffix = 'Complete', printEnd = "\r")

    pdf.close()
        
    # saving hits to json
    with open(output_path_json, 'w') as outjson:
        json.dump(best_hits, outjson)

    print(f"Total time: {time.time()-start:.2f}s")
    print(f"Average time to compute kmers: {np.mean(times_kmers):.2f}s")
    print(f"Average time to compute windows: {np.mean(times_windows):.2f}s")