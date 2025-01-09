from load_samples.ribo_db_build import prepare_ribo_db

import argparse

from ete3 import NCBITaxa
from taxoniq import Taxon

import sys, os, subprocess, shutil
from typing import Dict, List

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()

def getter(taxon : str, index : str, ribo_dir : str) -> None:
    """
    Run a subprocess to fetch the genome corresponding to the
    accession number.
    """
    subprocess.run(
            f'bio fetch {taxon} -t DNA -f fasta > "{index}/genome_{taxon}.fasta"',
            capture_output=True, shell=True
        )
    subprocess.run(
        f'bio fetch {taxon} -t DNA > "{ribo_dir}/ribosome_{taxon}.gb"',
            capture_output=True, shell=True
    )
    return 0


def parse_file(path : str) -> List[str]:
    """ Load the content of a file into a list of strings. """
    with open(path, 'r') as file:
        result = [line.strip('\n') for line in file.readlines()]
    return result

def taxname2taxid(taxlist : List[str]) -> Dict[str, List[str]]:
    """ Create a dictionnary from a list of species names """
    return NCBITaxa().get_name_translator(taxlist)

def taxid2seqid(taxid : str) -> str:
    """ Return the genome accession ID for a given taxid. """
    return Taxon(taxid).refseq_genome_accessions[0].accession_id

def generate_report(db_path : str, failed : list) -> None:
    """ Generates the downloading report. """
    with open(os.path.join(db_path, "_failed.txt"), 'w') as fail_file:
        for elt in failed:
            fail_file.write(f"{elt}\n")
    return 0

def main(taxon_file, db, ribo_dir, taxnames=True):
    samples = parse_file(taxon_file)
    if taxnames:
        taxon_dico = taxname2taxid(samples)   # {taxon name: taxon id}
    else:
        taxon_dico = {sample:sample for sample in samples} # {id: id}

    n = 0
    e = 0
    failed = []

    print("Starting Genomes Download...\n")
    progressbar(n, len(taxon_dico))

    for index, taxon in taxon_dico.items():
        n+=1
        enddir = os.path.join(db, f"{index.split('/')[0]}")
        os.makedirs(enddir)
        try :
            if taxnames:
                seq_id = taxid2seqid(taxon[0])
            else: 
                seq_id = taxon
            getter(seq_id, enddir, ribo_dir)
            progressbar(n, len(taxon_dico))

        except Exception:
            e+=1
            failed.append(index)
            os.rmdir(enddir)

    print(f"{n-e} taxa genomes download successful - {e} failed to download or unavailable\n")
    generate_report(db, failed)

    print("\n   Building Ribosomic sequence database\n")
    if not taxnames:
        n = prepare_ribo_db(ribo_dir, accession_id_flag = True)       
    else:
        n = prepare_ribo_db(ribo_dir, accession_id_flag = False)

    print(f"{n} Ribosomal informations gathered - {len(taxon_dico)-e-n} failed to parse.")
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download genomes from NCBI")
    parser.add_argument("taxon_file", help="File containing the list of taxons to download", type=str)
    parser.add_argument("-s", "--taxnames", help="Select if the list contains sequence ids and not taxon names", action="store_false")

    args = parser.parse_args()
    taxon_file = args.taxon_file
    taxon_basename = os.path.basename(taxon_file)
    taxnames = args.taxnames

    root_dir = os.path.abspath(os.path.join(__file__, "..", "..", ".."))

    if taxon_file == taxon_basename:
        taxon_file = os.path.join(root_dir, "input", "samples_list", taxon_basename)
    print(f"Taxon file: {taxon_file}")

    out_dir = os.path.join(root_dir, "input", "sequence_db", taxon_basename.split('.')[0])
    print(f"Output directory: {out_dir}")

    ribo_dir = os.path.join(root_dir, "input", "ribosome_db", taxon_basename.split('.')[0])
    print(ribo_dir)
    
    if os.path.exists(out_dir):
        print(f"Output directory already exists. If you continue, the content will be erased.")
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

    main(taxon_file, out_dir, ribo_dir, taxnames)

