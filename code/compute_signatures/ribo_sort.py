import os, re
from typing import List, Dict, Tuple, Generator
from collections import Counter

from compute_signatures.kmers import stream_kmers

from Bio import SeqIO


#################################################################################
#                         When the database is available                        #
#################################################################################

def parse_ribo_position(path_ribo_file : str) -> List[Tuple[int, int]]:
    """
    Parses the ribosomes database for their position in the sequence.
    """
    with open(path_ribo_file, 'r') as ribo_fasta:
        ribo_fasta = ribo_fasta.read()

    ribo_position = re.findall(r'>.*?(?<=\|).*?(?=\|).(.*?(?=\n))', ribo_fasta)
    ribo_position = [(int(elt.split('-')[0]), int(elt.split('-')[1])) for elt in ribo_position]

    return ribo_position

def drop_ribo_position(hits_position : Dict[str, float], 
                       path_ribo_db : str,
                       ) -> Dict[str, float]:
    """
    Filter out the positions corresponding to ribosome from the input hit dictionary, i.e. the dictionary 
    containing the position of incoherent signatures. 
    Do so by using the ribosome database computed when the genomes where downloaded. 
    Can be used only if the ribosomic data has been parsed previously. Otherwise, use the ribosomic signature module.
    """
    ribo_position_list = parse_ribo_position(path_ribo_db)

    position_map = list(map(lambda x : any(int(x) > y[0]-1000 and int(x) < y[1]+1000 for y in ribo_position_list), 
                    hits_position.keys())) #True if ribosome, False if not

    new_hits_position = {elt[0]:elt[1] for index, elt in enumerate(hits_position.items()) if position_map[index] == False}

    return new_hits_position

