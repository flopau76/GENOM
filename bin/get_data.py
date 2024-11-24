import get_taxid

import sys, os, subprocess

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r"):
    percent = f'{round(100 * (iteration / float(total)), 1)}'

    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)

    if iteration == total: 
        print()

def getter(taxon, index):
    os.chdir(index)
    test = subprocess.run(
            f"bio fetch {taxon} -t DNA -f fasta > genome_{taxon}.fasta",
            capture_output=True, shell=True
        )

    return 0

def load_taxid(taxon_file : str):
    return get_taxid.get_id(path=taxon_file)

def generate_report(db_path : str, failed : list):
    with open(os.path.join(db_path, "_failed.txt"), 'w') as fail_file:
        for elt in failed:
            fail_file.write(f"{elt}\n")

    return 0

def main(taxon_file, db):
    taxon_dico = load_taxid(taxon_file)
    n = 0
    e = 0
    failed = []
    for index, taxon in taxon_dico.items():
        enddir = os.path.join(db, f"{index}")
        os.makedirs(enddir)
        
        try :
            taxon = get_taxid.get_accession(taxon[0])
            getter(taxon, enddir)
            progressbar(n, len(taxon_dico))

            n+=1
        except Exception:
            n+=1
            e+=1
            failed.append(index)
            os.rmdir(enddir)
    
    generate_report(db, failed)
    print(f"{len(taxon_dico)-e} taxa genomes download successful\n{e} failed to download or unavailable")
    return 0

if __name__ == '__main__':
    taxon_file, db = sys.argv[1], sys.argv[2]
    
    if "report.txt" in os.listdir():
        os.remove(os.path.join(db, "report.txt"))

    main(taxon_file, db)

