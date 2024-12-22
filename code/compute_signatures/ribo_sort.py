import os, re
from typing import List, Dict, Tuple, Generator
from collections import Counter

import compute_signatures.ribo_class as rb
from compute_signatures.kmers import stream_kmers

from Bio import SeqIO


#################################################################################
#                         When the database is available                        #
#################################################################################

def parse_ribo_position(path_ribo_file : str) -> Generator[Dict[str, List[Tuple[int, int]]], any, any]:
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


################################################################################
#                       To remove if no use case is found                      #
################################################################################

def stream_ribo(ribo_db_path : str):
    """
    Passes through all ribosomes files and sort them into their corresponding 
    category.
    Separates Eukaryotes, Prokaryotes and Archea.
    """
    for file in os.listdir(ribo_db_path):
        file_path = os.path.join(ribo_db_path, file)
        all_ribo = [ribo_record for ribo_record in SeqIO.parse(file_path, "fasta")]

        yield all_ribo
            
def ribo_sorter(ribo_db_path : str) -> Dict:
    sender = {
        '5S': rb.ribo5S(kmer_count=Counter(), num=0),
        '16S': rb.ribo16S(kmer_count=Counter(), num=0),
        '23S': rb.ribo23S(kmer_count=Counter(), num=0),
        'SSU': rb.riboSSU(kmer_count=Counter(), num=0),
        'BSU': rb.riboBSU(kmer_count=Counter(), num=0)
        }
    
    for ribo_record_list in stream_ribo(ribo_db_path):
        for record in ribo_record_list:
            kmer_list = [kmer for kmer in stream_kmers([record.seq], k=8)]
            sender[record.id].num += len(kmer_list)
            sender[record.id].kmer_count = sender[record.id].kmer_count + Counter(kmer_list)

    print(sender)
    return


