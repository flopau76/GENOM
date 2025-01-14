import argparse

import sys, os, subprocess, shutil
from typing import Dict, List

from load_samples.get_taxid import parse_file
from load_samples.ribo_db_build import prepare_ribo_db

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()

def getter(seq_id : str, seq_dir : str, ribo_dir : str) -> None:
    """
    Run a subprocess to fetch the genome corresponding to the
    accession number.
    """
    subprocess.run(
            f'bio fetch {seq_id} -t DNA -f fasta > "{seq_dir}/{seq_id}.fasta"',
            capture_output=True, shell=True
        )
    subprocess.run(
        f'bio fetch {seq_id} -t DNA > "{ribo_dir}/{os.path.basename(seq_dir)}.gb"',
            capture_output=True, shell=True
        )
    return 0


def generate_report(db_path : str, failed : list) -> None:
    """ Generates the downloading report. """
    with open(os.path.join(db_path, "_failed.txt"), 'w') as fail_file:
        for elt in failed:
            fail_file.write(f"{elt}\n")
    return 0

def main(taxon_file, db, ribo_dir):
    n = 0
    e = 0

    taxon_dico, failed = parse_file(taxon_file)

    print("Starting Genomes Download...\n")
    progressbar(n, len(taxon_dico))

    for name, seq_id in taxon_dico.items():
        n+=1
        seq_dir = os.path.join(db, name)
        os.makedirs(seq_dir)
        try :
            getter(seq_id, seq_dir, ribo_dir)
        except Exception:
            e+=1
            failed.append(name)
            os.rmdir(seq_dir)
        progressbar(n, len(taxon_dico))

    print(f"{n-e} taxa genomes download successful - {len(failed)} failed to download or unavailable\n")
    generate_report(db, failed)

    print("\n   Building Ribosomic sequence database\n")
    n = prepare_ribo_db(ribo_dir)

    print(f"{n} Ribosomal informations gathered - {len(taxon_dico)-e-n} failed to parse.")
    return 0

if __name__ == '__main__':
    root_dir = os.path.abspath(os.path.join(__file__, "..", "..", ".."))
    parser = argparse.ArgumentParser(description="Download genomes from NCBI")
    parser.add_argument("taxon_file", help="File containing the list of taxons or sequence accesion ids to download", type=str)
    args = parser.parse_args()
    
    taxon_file = args.taxon_file
    taxon_basename = os.path.basename(taxon_file)
    if taxon_file == taxon_basename:
        taxon_file = os.path.join(root_dir, "input", "samples_list", taxon_basename)

    out_dir = os.path.join(root_dir, "input", "sequence_db", taxon_basename.split('.')[0])
    ribo_dir = os.path.join(root_dir, "input", "ribosome_db", taxon_basename.split('.')[0])

    # Check if the output directory already exists
    if os.path.exists(out_dir):
        print("Output directory already exists. If you continue, the content will be erased.")
        user_input = input("Do you want to continue? (Y/N): ")
        if user_input.strip().upper() != 'Y':
            print("Operation cancelled by user.")
            sys.exit(1)
        else:
            shutil.rmtree(out_dir)
            if os.path.exists(ribo_dir):
                shutil.rmtree(ribo_dir)

    os.makedirs(out_dir)
    os.makedirs(ribo_dir)

    main(taxon_file, out_dir, ribo_dir)
