import os, argparse
from typing import List
from time import time
import numpy as np

import analysis.loading as load
import analysis.origin_search as origin
import analysis.eval as eval
import analysis.graph as graph

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
    parser = argparse.ArgumentParser(description="Analyse the computed signatures and backtrack their origine from the database.")

    parser.add_argument("json_transfer_summary", help='The Path to the summary in .json format generated by the genome signature computation program. Must be in ./output/transfer_summary/')
    parser.add_argument("database_name", help="Path to the directory to analyze (whom signature have been computed). Must be in ./input/sequence_db/ or ./output/output_generator/")
    parser.add_argument("-o", '--output_file', help="Choose the name of the output file. Its destination is ./output/analysis/", type=str, default="analysis_report.txt")
    parser.add_argument('-r', '--reference_file', help="Path to the control file. Evaluate the results of the backtracking respect to it.", type=str, default=None)
    parser.add_argument('-w','--window_size', help="Choose the window size. Default length is 5000 (same as for signature computation).", type=int, default=10000)
    parser.add_argument('-n', '--number_kmer', help="Choose number of search kmers for backtracking. Default 3.", type=int, default=3)
    parser.add_argument('-t', '--threshold', help="Define the proportion of search kmer found in a genome to consider it as a hit. Default 60%", type = int, default=0.6)
    parser.add_argument('-p', '--probability', help='Define the maximum probability of finding a kmer of length n in a genome of size L at random. Default 0.1.', type=float, default=0.1)
    
    args = parser.parse_args()

    json_file = args.json_transfer_summary
    db_dir_name = args.database_name
    report_name = args.output_file
    eval_file_path = args.reference_file

    window_size = args.window_size
    number_search_kmer = args.number_kmer
    hit_threshold = args.threshold
    max_proba = args.probability

    time_ = []

    dir_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    json_path = os.path.join(dir_path, "output", "transfer_summary", json_file)
    db_path = os.path.join(dir_path, "input", "sequence_db", db_dir_name)
    output_path = os.path.join(dir_path, "output", "analysis")
    report_output_path = os.path.join(output_path, report_name)

    print("Loading transfer summary")
    ctrl_removal(output_path, report_name)

    init_report(report_output_path)

    print("Starting backtracking...")
    all_transfer = {}
    for transfer_summary in load.serialize_files(db_path, json_path=json_path):
        st = time()
        
        list_transfer = origin.find_kmer(db_path, transfer_summary, probability=max_proba, number_kmer=number_search_kmer, threshold=hit_threshold)
        write_report(report_output_path, list_transfer)

        time_.append(time()-st)

    print("Average runtime per transfer summary", round(np.average(time_)))
    print("Total Runtime", round(sum(time_)))
    
    if eval_file_path is not None:
        TP_rate, valid_list = eval.compare_files(report_output_path, eval_file_path, window_size=window_size)
    
        print(f"{round(TP_rate, 2)}% of the HGT from the reference file were found.")

        eval_report_name = os.path.join(output_path, f"eval_{report_name}")
        os.remove(eval_report_name)
        write_report(eval_report_name, valid_list)


"""
    print("Starting graph computation")
    out_graph = graph.dumy_dataset(report_output_path)


    plot_file = os.path.join(output_path, "graph_plot.png")
    graphml_file = os.path.join(output_path,"graph.graphml")
    degree_dist_file = os.path.join(output_path, "degree_distribution.png")
    graph.analyse_graph(out_graph, plot_file=plot_file, graphml_file=graphml_file, degree_dist_file=degree_dist_file)
"""