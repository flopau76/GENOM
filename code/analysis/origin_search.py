from Bio.Blast import NCBIWWW
from Bio import SeqIO
from analysis.loading import StrainHorizontalTransfer, HorizontalTransfer

from compute_signatures.kmers import stream_kmers, stream_kmers_file
import compute_signatures.metrics as metrics
from compute_signatures.display import display_windows
from compute_signatures.loading import open_genome

from scipy.special import kl_div
from typing import List, Dict, Tuple
from collections import Counter
from io import TextIOWrapper
import numpy as np
import os, time
from dataclasses import dataclass

@dataclass
class Conclusion:
    sender_found: str
    position_sender: str
    receiver : str
    position_receiver: str

def bootstrap_genome(seq : str, num_windows : int = 100, window_size : int = 2000, k_mer_length : int = 8) -> Dict[int, float]:
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

def load_bootstrapped_db(path_db : str, window_size : int = 2000, kmer_size : int = 8, Boostrap_iter : int = 100) -> Dict[str, Dict[int, float]]:
    """
    Load the boostrapped version of all the genomes in a given directory.
    Return their signature in a dictionnary
    """
    all_BS_db = {}
    for strain_dir in os.listdir(path_db):

        strain_dir_path = os.path.join(path_db, strain_dir)
        strain_file = os.listdir(strain_dir_path)[0]
        strain_path = os.path.join(strain_dir_path, strain_file)

        for strain_content in SeqIO.parse(strain_path, "fasta"):
            BSed_signature = bootstrap_genome(strain_content.seq, 
                                              num_windows=Boostrap_iter, 
                                              window_size=window_size, 
                                              k_mer_length=kmer_size)

        all_BS_db[strain_dir] = BSed_signature
    return all_BS_db

def screen_origins(list_transfer : StrainHorizontalTransfer, 
                   BSed_signature : Dict[str, Dict[int, float]],
                   kmer_size : int = 8) -> str:
    """
    Screens origin of a hit by computing the average signature over the whole genome
    for each strain in the database. Keeps closest one for window-analysis.
    """
    divergence_dico = {tr.start_position : (np.inf, 'NaN') for tr in list_transfer.transfer_summary}

    for strain_sender, dico_BS in BSed_signature.items():
        if strain_sender == list_transfer.strain:
            continue
        for transfer in list_transfer.transfer_summary:
            kmers_list = [kmer for kmer in stream_kmers([transfer.seq], kmer_size)]
            kmers_count = Counter(kmers_list)
            kmers_freq = {kmer: kmers_count[kmer]/len(kmers_list) for kmer in dico_BS.keys()}

            array_ref_Bootstrapped = np.array(list(dico_BS.values()))
            array_hit_sequence = np.array(list(kmers_freq.values()))
            
            kl = sum(kl_div(array_hit_sequence, array_ref_Bootstrapped)) #!! utilisation de la kl-div de scipy, à vérifier
            if kl < divergence_dico[transfer.start_position][0]:
                divergence_dico[transfer.start_position] = (kl, strain_sender)
        
    return divergence_dico

def KL_fixed_window_distance(fixed_window : List[int], file_pointer : TextIOWrapper, kmer_size = 8) -> np.ndarray:
    """
    Computes Kullback-Leibler divergence over the target genome of the possible transfer sequence with
    sequences of identical size. Idea is to find if there is a place in the genome with a similar profile.

    Query is the receiver's window .
    Target is the sender's genome.
    """
    # iterate once over the file to compute the total kmer frequency
    Target_list = list(stream_kmers_file(file_pointer, k=kmer_size))
    Target_count = set(Target_list)

    fixed_window_count = Counter(fixed_window)
    Query_fixed_profile = {kmer: (1+fixed_window_count[kmer])/(len(fixed_window)+len(Target_count)) for kmer in Target_count}

    distances = metrics.old_KL(kmers_list=Target_list, kmers_ref_freq=Query_fixed_profile, window_size=len(fixed_window))
    return distances


def sliding_window_search(top_screen : Dict[str, Dict[int,Tuple[float, str]]], 
                          db_path : str, 
                          transfer_summary : StrainHorizontalTransfer,
                          kmer_size : int = 8) -> List[Conclusion]:
    """
    Similarity measure between the transfered sequence found and a sliding window
    along the most likely strain to originate from after boostrapping.
    """
    target = top_screen[transfer_summary.strain]
    all_conclusion = []
    for index, transfer in enumerate(transfer_summary.transfer_summary):
        Query_sequence = [kmer for kmer in stream_kmers([transfer.seq], kmer_size)]
        Target_file = os.listdir(os.path.join(db_path, target[transfer.start_position][1]))[0]
        Target_pointer = os.path.join(db_path, target[transfer.start_position][1], Target_file)
        
        file_pointer = open_genome(Target_pointer)
        all_kl_div = KL_fixed_window_distance(Query_sequence, file_pointer)

        ccl = Conclusion(
            sender_found=target[transfer.start_position][1],
            position_sender=all_kl_div.argmin(),
            receiver= transfer_summary.strain,
            position_receiver= list(target.keys())[index]
        )

        all_conclusion.append(ccl)

    return all_conclusion
    

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


