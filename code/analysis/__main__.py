import os
from typing import List
from time import time
import numpy as np

import analysis.loading as load
import analysis.origin_search as origin

def init_report(path_report : str):
    """
    Initialize report file with fields.
    """
    with open(path_report, 'a+') as file:
        file.write(f"Sp_sending\t\tSending_position\t\tSp_receiving\t\tReceiving_position\n")

    return 0

def write_report(path_report : str, transfered : List[origin.Conclusion]):
    """
    Write the report following initialized format.
    """
    with open(path_report, 'a+') as file:
        for conclusion in transfered :
            file.write(f"{conclusion.sender_found}\t\t{conclusion.position_sender}\t\t{conclusion.receiver}\t\t{conclusion.position_receiver}\n")
    return 0

def ctrl_removal(base_dir : str, report_file : str):
    """
    Check for already present outputs and delete them to avoid error throwing.
    """
    if report_file in os.listdir(base_dir):
        os.remove(os.path.join(base_dir, report_file))
        print("  Output directory cleaned.\n")

    return 0

if __name__ == "__main__":
    
    json_file = "transfer_summary_data_test.json"
    db_dir_name = "data_test"
    report_name = "analysis_report.txt"
    time_ = []

    dir_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    json_path = os.path.join(dir_path, "output", json_file)
    db_path = os.path.join(dir_path, "input", "sequence_db", db_dir_name)
    output_path = os.path.join(dir_path, "output", "analysis")
    report_output_path = os.path.join(output_path, report_name)

    print("Loading transfer summary")
    ctrl_removal(output_path, report_name)
    init_report(report_output_path)

    print("Boostrapping Database Signatures...")
    Bootstrapped_signature = origin.load_bootstrapped_db(db_path, window_size=2000, kmer_size=8, Boostrap_iter=100)
    for transfer_summary in load.serialize_files(db_path, json_path=json_path):
        st = time()

        top_screen = {transfer_summary.strain : origin.screen_origins(transfer_summary, Bootstrapped_signature)}
        best_hit_search = origin.sliding_window_search(top_screen, db_path, transfer_summary)

        write_report(report_output_path, best_hit_search)
        time_.append(time()-st)
        break

    print("Average runtime per genome", round(np.average(time_)))
    print("Total Runtime", round(sum(time_)))
    