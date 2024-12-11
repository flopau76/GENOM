from Bio.Blast import NCBIWWW
from Bio import SeqIO
from analysis.loading import StrainHorizontalTransfer, HorizontalTransfer

from compute_signatures.kmers import stream_kmers
from compute_signatures.signatures import KLdivergence
from compute_signatures.display import display_windows

from scipy.special import kl_div
from typing import List, Dict, Tuple
from collections import Counter
import numpy as np
import os

def bootstrap_genome(seq : str, num_windows : int = 100, window_size : int = 2000, k_mer_length : int = 8):
    """
    Randomly select n overlapping-able windows of length w and computes
    an average signature of k-mer frequency of length k.
    By default the parameters are the following : 
        n = 100
        w = 2000
        k = 8
    """
    start_array = np.random.randint(0, len(seq), size=num_windows)
    seq_list = [seq[i:i+window_size] for i in start_array]
    kmers_list = [kmer for kmer in stream_kmers(seq_list, k_mer_length)]
    
    BS_kmer_freq = {kmer: count/(len(kmers_list)) for kmer, count in Counter(kmers_list).items()}
    return BS_kmer_freq

def screen_origins(path_db : str, 
                   list_transfer : StrainHorizontalTransfer, 
                   window_size : int = 2000, 
                   kmer_size : int = 8) -> str:
    """
    Screens origin of a hit by computing the average signature over the whole genome
    for each strain in the database. Keeps closest one for window-analysis.
    """
    divergence_dico = {tr.start_position : (np.inf, 'NaN') for tr in list_transfer.transfer_summary}
    for strain_dir in os.listdir(path_db):
        
        if strain_dir == list_transfer.strain or strain_dir.endswith('.txt') or strain_dir == 'ribo_db':
            continue

        strain_dir_path = os.path.join(path_db, strain_dir)
        strain_file = os.listdir(strain_dir_path)[0]
        strain_path = os.path.join(strain_dir_path, strain_file)

        for strain_content in SeqIO.parse(strain_path, "fasta"):
            BSed_signature = bootstrap_genome(strain_content.seq, window_size=window_size, k_mer_length=kmer_size)

        for transfer in list_transfer.transfer_summary:
            kmers_list = [kmer for kmer in stream_kmers([transfer.seq], kmer_size)]
            kmers_count = Counter(kmers_list)
            kmers_freq = {kmer: kmers_count[kmer]/len(kmers_list) for kmer in BSed_signature.keys()}

            array_ref_Bootstrapped = np.array(list(BSed_signature.values()))
            array_hit_sequence = np.array(list(kmers_freq.values()))
            
            kl = sum(kl_div(array_hit_sequence, array_ref_Bootstrapped))
            if kl < divergence_dico[transfer.start_position][0]:
                divergence_dico[transfer.start_position] = (kl, strain_dir)
        
    return divergence_dico

def load_target_sequence(dir_name : str, db_path : str, kmer_size : int = 8) -> List[str]:
    """
    Load as a list of kmer of size k the possible origin of the transfer previously found.
    Return it in as the reference to be the used for the sliding window comparison.
    """
    target_dir_path = os.path.join(db_path, dir_name)
    target_file_path = os.path.join(target_dir_path, os.listdir(target_dir_path)[0])

    for content in SeqIO.parse(target_file_path, "fasta"):
        seq = content.seq
    return [kmer for kmer in stream_kmers([seq], kmer_size)]

def KL_fixed_window_distance(fixed_window : List[int], Target_kmers_list : List[int]) -> np.ndarray:
    """
    Computes Kullback-Leibler divergence over the target genome of the possible transfer sequence with
    sequences of identical size. Idea is to find if there is a place in the genome with a similar profile.
    """
    fixed_window_count = Counter(fixed_window)
    Target_count = Counter(Target_kmers_list)
    Target_profile = {kmer : (1+count)/(len(Target_kmers_list)+len(Target_count)) for kmer, count in Target_count.items()}
    Query_fixed_profile = {kmer: (1+fixed_window_count[kmer])/(len(fixed_window)+len(Target_count)) for kmer in Target_profile.keys()}

    return KLdivergence(Target_kmers_list, Query_fixed_profile, window_size=len(fixed_window))


def sliding_window_search(top_screen : Dict[str, Dict[int,Tuple[float, str]]], 
                          db_path : str, 
                          transfer_summary : StrainHorizontalTransfer,
                          kmer_size : int = 8):
    """
    Similarity measure between the transfered sequence found and a sliding window
    along the most likely strain to originate from after boostrapping.
    """
    target = top_screen[transfer_summary.strain]
    for transfer in transfer_summary.transfer_summary:
        Query_sequence = [kmer for kmer in stream_kmers([transfer.seq], kmer_size)]
        Target_sequence = load_target_sequence(target[transfer.start_position][1], db_path, kmer_size)

        all_kl_div = KL_fixed_window_distance(Query_sequence, Target_sequence)
        display_windows(all_kl_div, f"KL divergence between a sequence from {transfer_summary.strain} at position {transfer.start_position}\nand the full genome of {target[transfer.start_position][1]} whom may be the origin of a horizontal transfer")
        return



def blast_seq(transfer_elt : HorizontalTransfer) -> str:
    """
    Removed for now.
    """
    result_handle = NCBIWWW.qblast("blastn", "nt", transfer_elt.seq, format_type="Text")
    text = result_handle.read()#.split("ALIGNMENTS")[0]
    return text

def merge_seq(list_transfer: List[HorizontalTransfer], window_size : int = 2000) -> List[HorizontalTransfer]:
    """
    Deprecated. Removed from __main__.
    """
    list_transfer = sorted(list_transfer, key= lambda x : x.start_position)
    cluster_transfer = []
    i = 0
    while i < len(list_transfer):
        cluster = [elt for elt in list_transfer 
                   if elt.start_position >= list_transfer[i].start_position 
                        and elt.start_position <= list_transfer[i].start_position + window_size]
        cluster_transfer.append(cluster)
        i += len(cluster)
        
    merged_list = []
    while cluster_transfer:
        cluster_elts = cluster_transfer.pop(0)
        start = cluster_elts[0].start_position
        end_start = cluster_elts[-1].start_position-1

        start_seq = cluster_elts[0].seq[0:end_start-start]
        merged_seq = HorizontalTransfer(
            ID=cluster_elts[0].ID,
            description=cluster_elts[0].description,
            start_position=start,
            end_position=cluster_elts[-1].end_position,
            divergence=cluster_elts[0].divergence,
            seq = start_seq + cluster_elts[-1].seq,
        )
        merged_list.append(merged_seq)

    return merged_list


