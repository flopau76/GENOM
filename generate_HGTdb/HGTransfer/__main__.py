import os, shutil
import numpy as np

from HGTransfer.load_input import HGT, loader
from HGTransfer.exchanger import transferer, write_output_db

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

def init_report(path_report : str):
    """
    Initialize report file with fields.
    """
    with open(path_report, 'a+') as file:
        file.write(f"Sp_sending\t\tSending_position\t\tSp_receiving\t\tReceiving_position\n")

    return 0

def write_report(path_report : str, transfered : HGT):
    """
    Write the report following initialized format.
    """
    sender = transfered.sender_object
    receiver = transfered.receiver_object
    with open(path_report, 'a+') as report:
        report.write(f"{sender.strain}\t\t{sender.transfer_start}-{sender.transfer_end}\t\t{receiver.strain}\t\t{receiver.reception_position}\n")

    return 0

def ctrl_removal(base_dir : str, out_dir : str, report_file : str):
    """
    Check for already present outputs and delete them to avoid error throwing.
    """
    path_outdir = os.path.join(base_dir, out_dir)
    if out_dir in os.listdir(base_dir):
        shutil.rmtree(path_outdir)
        os.makedirs(path_outdir)

    if report_file in os.listdir(base_dir):
        os.remove(os.path.join(base_dir, report_file))

    print("  Output directory cleaned.\n")
    return 0

if __name__ == "__main__":
    print("Starting HGTdb generation...")
    iterations = 1000
    transfer_proba = 0.05

    base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    in_dir, out_dir, report_file = "input", 'output', "HGT_report.txt"

    db_path = os.path.join(base_dir, in_dir)
    output_path = os.path.join(base_dir, out_dir)
    path_report = os.path.join(base_dir, report_file)

    if out_dir in os.listdir(base_dir) or report_file in os.listdir(base_dir):
        print("  Cleaning output directory...")
        ctrl_removal(base_dir, out_dir, report_file)

    init_report(path_report)
    for iteration in range(0, iterations):
        progressbar(iteration=iteration+1, total=iterations)
        test = np.random.random(1)
        if test < transfer_proba:
            selected = loader(db_path, iteration)
            transfered = transferer(selected, iteration, iterations)

            write_report(path_report, transfered)
            write_output_db(output_path, transfered)

    print(f"\nHGT database is Ready. Refer to the file {report_file} for the list of the transfers that occured.\n")
    

