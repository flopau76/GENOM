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
    with open(path_report, 'a') as file:
        file.write("Sp_sending\t\tSending_start\t\tSending_end\t\tSp_receiving\t\tReceiving_position\n")
    return 0

def write_report(path_report : str, transfered : HGT):
    """ Write the report following initialized format. """
    sender = transfered.sender_object
    receiver = transfered.receiver_object
    with open(path_report, 'a+') as report:
        report.write(f"{transfered.iteration} {sender.strain}\t\t{sender.transfer_start}\t\t{sender.transfer_end}\t\t{transfered.iteration} {receiver.strain}\t\t{receiver.reception_position}\n")

    return 0


if __name__ == "__main__":
    print("Starting HGTdb generation...")
    
    parser = argparse.ArgumentParser(description="Generate Horizontal Transfer from a provided database")
    parser.add_argument('-db', '--input_db', help="Input database with list of taxa", type=str, default ="db")
    parser.add_argument('-o', '--output_db', help="Name of the directory in `input/sequence_db` where generated data is stored in", type=str, default="generator_db")
    parser.add_argument('-r', '--report', help='Name of the report file. Is found in input/sequence_db/<output_db>', type=str,default='HGT_report.txt')
    parser.add_argument('-i', '--iterations', help="Number of iterations of transfer trials", type=int, default=1000)
    parser.add_argument('-p', '--probability', help="Set horizontal transfer probability", type=float, default=0.01)

    args = parser.parse_args()
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    in_dir = args.input_db
    out_dir = args.output_db
    report_file = args.report
    iterations = args.iterations
    transfer_proba = args.probability

    base_db = os.path.join(base_dir, "input", "sequence_db")
    input_path = os.path.join(base_db, in_dir) #send to input generator
    output_path_db = os.path.join(base_db, out_dir) #send to input as modified database 
    
    path_report = os.path.join(output_path_db, report_file)

    if os.path.exists(output_path_db):
        print("  Cleaning output directory...")
        shutil.rmtree(output_path_db)
        os.makedirs(output_path_db)
    
    taken = []
    init_report(path_report)
    for iteration in range(0, iterations):
        progressbar(iteration=iteration+1, total=iterations)
        test = np.random.random(1)
        if test < transfer_proba:
            selected, new_taken = loader(input_path, iteration, taken)
            transfered = transferer(selected, iteration, iterations)

            write_report(path_report, transfered)
            write_output_db(output_path_db, transfered, iteration)

            taken += new_taken
            taken = list(set(taken)) #2 files cannot be taken twice for HGT
            
    
    print(f"\nHGT database is Ready.\nRefer to the file {report_file} for the list of the transfers that occured.\nFiles can be found in the input/sequence_db/{out_dir}")
    

