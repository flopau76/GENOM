import load_samples.get_taxid as get_taxid
from load_samples.ribo_db_build import prepare_ribo_db

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
        f'bio fetch {taxon} -t DNA > "{ribo_dir}/ribosome_{taxon}.gb"',
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

def main(taxon_file, db, ribo_dir):
    taxon_dico = load_taxid(taxon_file)
    os.makedirs(out_dir)
    os.makedirs(ribo_dir)
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
            taxon = get_taxid.get_accession(taxon[0])
            getter(taxon, enddir, ribo_dir)
            progressbar(n, len(taxon_dico))
        except Exception:
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
    if len(sys.argv) != 2:
        print("Usage: python -m load_sample  <sample_file>")
        sys.exit(1)

    root_dir = os.path.abspath(os.path.join(__file__, "..", "..", ".."))

    taxon_file = sys.argv[1]
    taxon_basename = os.path.basename(taxon_file)
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

    main(taxon_file, out_dir, ribo_dir)

