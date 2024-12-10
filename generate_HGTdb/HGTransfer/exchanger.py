import os, random
import numpy as np
from typing import Union

from HGTransfer.load_input import HGT, Sender, Receiver
import HGTransfer.events as events

def transferer(selected : HGT, current_iter : int, total_iter : int) -> HGT:
    """
    Transfers a randomly selected sequence of consequent size from the sender 
    to the receiver organism.
    """
    sender = selected.sender_object
    receiver = selected.receiver_object

    seq_transfered = sender.seq[sender.transfer_start:sender.transfer_end]
    receiver.seq = receiver.seq[0:receiver.reception_position] + seq_transfered + receiver.seq[receiver.reception_position:len(receiver.seq)]

    transfered = HGT(
        sender_object=sender,
        receiver_object=receiver,
        iteration=selected.iteration
    )
    
    #transfered = add_noise(transfered, current_iter, total_iter)

    return transfered

def add_noise(transfered : HGT, current_iter : int, total_iter : int) -> HGT:
    """
    Adds some noise to the transfered sequences for realism (and challenge).
    The statistical chances of every event occuring were taken from the RecombSimulator
    program developed by Raphaël CHAMPEIMONT during his PhD and which we studied in PHYG.
    We consider an average of 2 events per iteration, with a standard deviation of k/5 with
    k the number of event left to occur.

    Citation :
    Raphael Champeimont. Genetic mutations combinatorics. Bioinformatics [q-bio.QM]. 
    Université Pierre et Marie Curie - Paris VI, 2014. English. ⟨NNT : 2014PA066636⟩. ⟨tel-01118660v2⟩
    """
    k = total_iter-current_iter
    event_number = round(np.random.normal(k*2, k/5))

    sender = transfered.sender_object
    receiver = transfered.receiver_object
    while event_number != 0:
        
        prop_SNP_sender, prop_SNP_receiver = np.random.random(1), np.random.random(1) 
        sender, receiver = events.SNPs(sender, prop_SNP_sender), events.SNPs(receiver, prop_SNP_receiver)

        event_number -=1
    return transfered

def noise_dispach(victim : Union[Sender, Receiver], roll : np.ndarray):
    for proba in roll:
        pass

def write_output_db(output_path_receiver : str, output_path_sender : str, transfer : HGT):
    """
    Writes the output genomes from the transfer in a new directory for further analysis 
    if required.
    Each output genome is written in a newly created directory named after the iteration in
    which the transfer happened. Files are named with the format <role>_<strain_name> and in fasta format.
    """
    sender = transfer.sender_object
    receiver = transfer.receiver_object

    output_dir_sender = os.path.join(output_path_sender, str(transfer.iteration))
    os.makedirs(output_dir_sender)
    output_dir_receiver = os.path.join(output_path_receiver, str(transfer.iteration))
    os.makedirs(output_dir_receiver)

    sender_path = os.path.join(output_dir_sender, f"sender_{sender.strain}.fasta")
    receiver_path = os.path.join(output_dir_receiver, f"receiver_{receiver.strain}.fasta")
    
    with open(sender_path, "w") as sender_file:
        sender_file.write(f">{sender.strain}_{sender.ID} | {sender.transfer_start}-{sender.transfer_end}\n{sender.seq}\n")
    
    with open(receiver_path, 'w') as receiver_file:
        receiver_file.write(f">{receiver.strain}_{receiver.ID} | {receiver.reception_position}\n{receiver.seq}\n")




