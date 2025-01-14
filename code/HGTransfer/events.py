from typing import Union
import numpy as np
import re

from HGTransfer.load_input import Sender, Receiver

def inverstion(victim : Union[Sender, Receiver]) -> Union[Sender, Receiver]:
    """
    Reverses random parts of the genome in either the sender or the receiver.
    """
    pass

def translocation(victim : Union[Sender, Receiver])-> Union[Sender, Receiver]:
    """
    Translocates and exchange random parts of the genome in either the sender or the receiver.
    """
    pass

def deletion(victim : Union[Sender, Receiver])-> Union[Sender, Receiver]:
    """
    Deletes random parts of the genome in either the sender or the receiver.
    """
    pass

def duplication(victim : Union[Sender, Receiver]) -> Union[Sender, Receiver]:
    """
    Duplicates random parts of the genome in either the sender or the receiver.
    """
    pass

def SNPs(victim : Union[Sender, Receiver], prop : np.ndarray, number: int = 100) -> Union[Sender, Receiver]:
    """
    Creates SNPs randomly in the genome of either the send or the receiver.
    We consider the number of SNPs to occur as prop*number where prop is the 
    fifth roll casted previously.

    Supposing equirepartition of SNPs types.
    """
    num_SNPs = int(round(prop[0]*number)/3)
    SNPs_types = {"insertion": insertion, "deletion":deletion, "substitution":substitution}

    for event, func in SNPs_types.items():
        victim.seq = func(victim.seq, num_SNPs)
    
    return victim

def insertion(seq : str, num : int) -> str:
    """
    Insert num elements in the sequence seq
    """
    seq = str(seq)
    insert_list = np.random.choice(["A", 'T', 'C', 'G'], size = num, replace=True)
    insert_positions = np.random.randint(low=0, high=len(seq), size=num)

    for index, position in enumerate(insert_positions):
        pattern = re.escape(seq[position])
        seq = re.sub(pattern, f"{seq[position]}{insert_list[index]}", seq, count=1)

    return seq

def deletion(seq : str, num : int) -> str:
    """
    Remove num elements from the sequence seq
    """
    seq = str(seq)
    del_positions = np.random.randint(low=0, high=len(seq), size=num)

    for position in del_positions:
        pattern = re.escape(seq[position])
        seq = re.sub(pattern, "", str(seq), count=1)

    return seq

def substitution(seq : str, num : int) -> str:
    """
    Substitute num elements from the sequence seq
    """
    seq = str(seq)
    insert_list = np.random.choice(["A", 'T', 'C', 'G'], size = num, replace=True)
    insert_positions = np.random.randint(low=0, high=len(seq), size=num)

    for index, position in enumerate(insert_positions):
        pattern = re.escape(seq[position])
        seq = re.sub(pattern, insert_list[index], str(seq), count=1)
    return seq