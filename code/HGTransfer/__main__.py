import os, shutil, argparse
import numpy as np

from HGTransfer.load_input import HGT, loader
from HGTransfer.exchanger import transferer, write_output_db


def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'
    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    if iteration == total: 
        print()

def init_report(path_report : str):
    """ Initialize report file with fields. """
    with open(path_report, 'a+') as file:
        file.write(f"Sp_sending\t\tSending_position\t\tSp_receiving\t\tReceiving_position\n")
    return 0

def write_report(path_report : str, transfered : HGT):
    """ Write the report following initialized format. """
    sender = transfered.sender_object
    receiver = transfered.receiver_object
    with open(path_report, 'a+') as report:
        report.write(f"{transfered.iteration} {sender.strain}\t\t{sender.transfer_start}-{sender.transfer_end}\t\t{transfered.iteration} {receiver.strain}\t\t{receiver.reception_position}\n")

    return 0


if __name__ == "__main__":
    print("Starting HGTdb generation...")
    iterations = 1000
    transfer_proba = 0.01

    parser = argparse.ArgumentParser(description="Generate Horizontal Transfer from a provided database")
    parser.add_argument('input_db', help="Name of the database containing fasta sequences (must be in input/sequence_db)", type=str)
    args = parser.parse_args()
    in_dir = args.input_db

    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    input_path = os.path.join(base_dir, "input", "sequence_db", in_dir)
    output_path = os.path.join(base_dir, "input", "generated_HGT", in_dir)

    output_path_receiver = os.path.join(output_path, "output_receiver")
    output_path_sender = os.path.join(output_path, "output_sender")
    output_path_report = os.path.join(output_path, "HGT_report.txt")

    if os.path.exists(output_path):
        print("  Cleaning output directory...")
        shutil.rmtree(output_path)
        os.makedirs(output_path)
    os.makedirs(output_path_receiver)
    os.makedirs(output_path_sender)

    init_report(output_path_report)
    for iteration in range(0, iterations):
        progressbar(iteration=iteration+1, total=iterations)
        test = np.random.random(1)
        if test < transfer_proba:
            selected = loader(input_path, iteration)
            transfered = transferer(selected, iteration, iterations)

            write_report(output_path_report, transfered)
            write_output_db(output_path_receiver, output_path_sender, transfered, iteration)

    print(f"\nHGT database is Ready. Refer to the file {output_path_report} for the list of the transfers that occured.\n")
    

