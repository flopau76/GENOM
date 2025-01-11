import json, os
from typing import Dict, Generator
from dataclasses import dataclass

from Bio import SeqIO


@dataclass
class StrainHorizontalTransfer:
    strain : str
    transfer_summary : list

@dataclass
class HorizontalTransfer:
    ID : str
    description : str
    start_position : int
    end_position : int
    divergence : float
    seq : str

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """
    Show a progress bar indicating downloading progress
    """
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()

def load_json(path : str) -> Dict[str, Dict[str, float]]:
    with open(path, 'r', encoding='utf-8') as summary:
        dico = json.load(summary)
    return dico

def serialize_files(dir_path : str,
                    json_path : str = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"transfer_summary.json"),
                    window_size : int = 5000) -> Generator:

    strain_dico = load_json(json_path)

    n = 0
    for strain, position_dico in strain_dico.items():
        n+=1
        progressbar(n, len(strain_dico))
        directory = os.path.join(dir_path, strain)
        file = os.listdir(directory)[0]
        for seq in SeqIO.parse(os.path.join(directory, file), 'fasta'):
            list_transfer = []
            for position, divergence in position_dico.items():
                position = int(position)
                window = seq.seq[position:position+window_size]
                list_transfer.append(HorizontalTransfer(seq.id, 
                                                        seq.description,
                                                        position, 
                                                        position+window_size,
                                                        divergence,
                                                        seq=window))
        
        yield StrainHorizontalTransfer(
            strain,
            transfer_summary=list_transfer
        )
        
                

