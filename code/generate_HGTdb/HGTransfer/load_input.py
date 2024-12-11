import os
import numpy as np
from dataclasses import dataclass

from Bio import SeqIO

@dataclass
class Receiver:
    strain : str
    strain_file : str
    ID : str
    seq : str
    origin : str
    reception_position : int

@dataclass
class Sender:
    strain : str
    strain_file : str
    ID : str
    seq : str
    destination: str
    transfer_start : int
    transfer_end : int

@dataclass
class HGT:
    sender_object : Sender
    receiver_object : Receiver
    iteration : int


def selector(input_path : str) -> np.ndarray[str]:
    """
    Randomly select two directory in the input directory where all strain directories are kept. 
    """
    return np.random.choice(os.listdir(input_path), size=2)

def loader(input_path : str, iteration : int) -> HGT:
    """
    Load both files previously selected. One as the sender and the other as the receiver.
    Return an HGT (Horizontal Gene Transfer) object containing both Sender and Receiver objects.

    Also randomly select the boundaries of the sequence transfered from the sender and the position 
    to which the sequence is transfered to the receiver.
    """
    files = selector(input_path)
    sender = files[0]
    receiver = files[1]

    sender_dir = os.path.join(input_path, sender)
    receiver_dir = os.path.join(input_path, receiver)

    sender_file = os.path.join(sender_dir, os.listdir(sender_dir)[0])
    receiver_file = os.path.join(receiver_dir, os.listdir(receiver_dir)[0])

    for record in SeqIO.parse(sender_file, "fasta"):
        length_transfer = np.random.randint(low=1000, high=20000)
        transfer_start=np.random.randint(0, len(record.seq)-length_transfer)
        transfer_end = transfer_start+length_transfer

        sender_obj = Sender(
            strain=sender,
            strain_file=sender_file,
            ID=record.id,
            seq=record.seq,
            destination=receiver,
            transfer_start=transfer_start,
            transfer_end=transfer_end
        )

    for record in SeqIO.parse(receiver_file, "fasta"):
        receiver_obj = Receiver(
            strain=receiver,
            strain_file=receiver_file,
            ID=record.id,
            seq=record.seq,
            origin=sender,
            reception_position=np.random.randint(low=0, high=len(record.seq))
        )
    
    return HGT(
        sender_object=sender_obj,
        receiver_object=receiver_obj,
        iteration=iteration
    )