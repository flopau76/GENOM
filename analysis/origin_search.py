from Bio.Blast import NCBIWWW
from analysis.loading import HorizontalTransfer

from typing import List


def blast_seq(transfer_elt : HorizontalTransfer) -> str:
    result_handle = NCBIWWW.qblast("blastn", "nt", transfer_elt.seq, format_type="Text")
    text = result_handle.read().split("ALIGNMENTS")[0]
    return text

def merge_seq(list_transfer: List[HorizontalTransfer], window_size : int = 2000) -> List[HorizontalTransfer]:
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


