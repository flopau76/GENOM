import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
from scipy.signal import find_peaks

from collections import Counter
import argparse
import os, json, time

import compute_signatures.metrics as metrics
from compute_signatures.loading import iter_directory, open_genome
from compute_signatures.kmers import stream_kmers_file
import compute_signatures.display as display

from typing import Dict

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()

def find_potential_HGT(result:np.ndarray) -> Dict[int, float]:
    """ Find potential HGT regions corresponding to the highest peaks in the result """
    highest_values_indices, _ = find_peaks(result, prominence=(np.max(result) - np.min(result))/3)
    highest_values = result[highest_values_indices]
    best_hits = {int(idx):val for idx, val in zip(highest_values_indices, highest_values)}
    return best_hits

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the signature of a genome and find potential HGT regions")
    parser.add_argument("input_db", help="The name of the input database (must be in `input/sequence_db/`)")
    parser.add_argument('-k', '--kmer', help='The size of the kmer (default=8)', type=int, default=8)
    parser.add_argument('-w', '--window', help='The size of the sliding window (default=2000)', type=int, default=2000)

    args = parser.parse_args()

    k = args.kmer
    window_size = args.window
    folder_name = args.input_db

    base_dir =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    input_folder = os.path.join(base_dir, "input", "sequence_db", folder_name)
    output_path = os.path.join(base_dir, "output", f"transfer_summary_{folder_name}.json")
    output_path_pdf = os.path.join(base_dir, "output", f"transfer_summary_{folder_name}.pdf")
    tmp_dir_images = os.path.join(base_dir, "output", f"images_{folder_name}")

    if not os.path.exists(tmp_dir_images):
        os.makedirs(tmp_dir_images)

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
        metric = metrics.Convolution(kmers_freq)
        result = metric.slide_window(kmers_list, window_size)
        t2 = time.time()
        times_windows.append(t2-t1)

        # find the highest peaks
        sample_hits = find_potential_HGT(result)
        best_hits[sample] = sample_hits

        # save the resulting figure
        fig = display.display_windows(result, hits=sample_hits, title=f"{sample}", ylabel="KL divergence", dpi=300)
        fig.savefig(os.path.join(tmp_dir_images, f"{sample}.png"))
        plt.close("all")

        n += 1
        progressbar(n, n_total , prefix = 'Progress:', suffix = 'Complete', printEnd = "\r")

    # saving hits to json
    with open(output_path, 'w') as outjson:
        json.dump(best_hits, outjson)

    # Combine all the images into a single PDF
    with bpdf.PdfPages(output_path_pdf) as pdf:
        for sample in best_hits.keys():
            fig_path = os.path.join(tmp_dir_images, f"{sample}.png")
            if os.path.exists(fig_path):
                fig = plt.figure()
                img = plt.imread(fig_path)
                plt.imshow(img)
                plt.axis('off')
                pdf.savefig(fig)
                plt.close(fig)
                os.remove(fig_path)
    os.rmdir(tmp_dir_images)

    print(f"Total time: {time.time()-start:.2f}s")
    print(f"Average time to compute kmers: {np.mean(times_kmers):.2f}s")
    print(f"Average time to compute windows: {np.mean(times_windows):.2f}s")