import get_taxid
from ribo_db_build import prepare_ribo_db

import sys, os, subprocess, shutil
from typing import Dict

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

def getter(taxon : str, index : str, ribo_dir : str) -> None:
    """
    Run a subprocess to fetch the genome corresponding to the
    accession number.
    """
    test = subprocess.run(
            f'bio fetch {taxon} -t DNA -f fasta > "{index}/genome_{taxon}.fasta"',
            capture_output=True, shell=True
        )
    subprocess.run(
        f'bio fetch {taxon} -t DNA > "{ribo_dir}/ribosome_{taxon}.gb',
            capture_output=True, shell=True
    )
    return 0

def load_taxid(taxon_file : str) -> Dict[str, str]:
    """
    Load dictionnary associating a species name with its taxid
    """
    return get_taxid.get_id(path=taxon_file)

def generate_report(db_path : str, failed : list) -> None:
    """
    Generates the downloading report.
    """
    with open(os.path.join(db_path, "_failed.txt"), 'w') as fail_file:
        for elt in failed:
            fail_file.write(f"{elt}\n")

    return 0

def main(taxon_file, db):
    taxon_dico = load_taxid(taxon_file)
    n = 1
    e = 0
    failed = []
    ribo_dir = os.path.join(db, 'ribo_db')

    os.makedirs(ribo_dir)

    print("Starting Genomes Download...\n")
    for index, taxon in taxon_dico.items():

        enddir = os.path.join(db, f"{index.split('/')[0]}")
            
        os.makedirs(enddir)

        try :
            taxon = get_taxid.get_accession(taxon[0])
            getter(taxon, enddir, ribo_dir)
            progressbar(n, len(taxon_dico))

            n+=1
        except Exception:
            n+=1
            e+=1
            failed.append(index)
            os.rmdir(enddir)
    
    print("\n   Building Ribosomic sequence database\n")
    n = prepare_ribo_db(ribo_dir)
    generate_report(db, failed)

    print(f"{len(taxon_dico)-e} taxa genomes download successful - {e} failed to download or unavailable\n")
    print(f"{n} Ribosomal informations gathered - {len(taxon_dico)-e-n} failed to parse.")
    return 0

if __name__ == '__main__':
    taxon_file, db = sys.argv[1], sys.argv[2]
    
    if "report.txt" in os.listdir():
        os.remove(os.path.join(db, "report.txt"))

    main(taxon_file, db)

