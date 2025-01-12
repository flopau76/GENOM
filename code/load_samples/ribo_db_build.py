import re, os
from dataclasses import dataclass
from typing import List, Dict, Tuple, Generator

def progressbar(iteration, total, prefix = '', suffix = '', filler = '█', printEnd = "\r") -> None:
    """ Show a progress bar indicating downloading progress """
    percent = f'{round(100 * (iteration / float(total)), 1)}'
    add = int(100 * iteration // total)
    bar = filler * add + '-' * (100 - add)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    if iteration == total: 
        print()

@dataclass
class GenBankElement:
    organism_name : str
    organism_realm : str
    accession_id : str
    ribosome_beacon : Dict[str, Tuple[int, int]]
    ribosome_sequences : Dict[str, str]

def postprocess_regex(pattern : List[Tuple[str]]) -> List[List[str]]:
    liste = []
    for elt in pattern:
        new_tup = [i for i in elt if i != '']
        liste.append(new_tup)
    return liste

def parser(ribo_file_content : str):
    """
    Parses the content of a .gb using RegEx and return the corresponding GenBankElement.
    Note : sensible to the GenBank Summary format.
    """
    genome_seq = ''.join(re.findall('[atcg]+', ribo_file_content))
    accession_id = re.findall(r'VERSION     (.*?(?=\n))', ribo_file_content)[0]
    organism_name = re.findall(r'ORGANISM..(.+?)(?=\n)', ribo_file_content)[0]
    realm = re.findall(r'ORGANISM..[\s\S]*?(?<=\s{13})(.*?)(?=;)', ribo_file_content)[0]
        
    pattern = re.findall(r'rRNA\s+(?:complement\(join\(([\d\.\.,\s]+)\)\)|complement\((\d+\.\.\d+)\)|(\d+\.\.\d+))\s+.*\s+.*\s+.*?(?<=product=\")(.+?)(?=\")',
                    ribo_file_content)
    pattern = postprocess_regex(pattern)
    ribo_dico = {elt[1] : [] for elt in pattern}

    for elt in pattern:
        beacons, ribo_id = elt
        beacons_list = re.findall('([0-9]+..[0-9]+)', beacons)

        ribo_dico[ribo_id].extend([(int(beacons.split('..')[0]), int(beacons.split('..')[1])) for beacons in beacons_list])

    ribo_seq_dico = {ribo_id : [genome_seq[beacon[0]:beacon[1]] for beacon in beacons] for ribo_id, beacons in ribo_dico.items()}      
    
    return GenBankElement(
            ribosome_sequences=ribo_seq_dico,
            organism_realm=realm,
            organism_name=organism_name,
            ribosome_beacon=ribo_dico,
            accession_id=accession_id
        )

def save_to_fasta(ribo_elts : GenBankElement, out_path : str) -> None:
    """ Write ribosomic sequences to a new fasta file. """
    with open(out_path, 'a+') as ribofile:
        for ribo_id, beacons in ribo_elts.ribosome_beacon.items():
            for index, seq in enumerate(ribo_elts.ribosome_sequences[ribo_id]):
                ribofile.write(
                    f">{ribo_id}_{index}_{ribo_elts.organism_name}|{ribo_elts.organism_realm}|{beacons[index][0]}-{beacons[index][1]}\n{seq.upper()}\n"
                )

def prepare_ribo_db(ribo_db_dir : str):
    """ Transform all gb files of the directory into fasta file containing ribosomic sequence. """
    n = 0
    files = os.listdir(ribo_db_dir)
    progressbar(n, len(files))
    for file in files:
        n+=1
        if not file.endswith('.gb'):
            continue
        path = os.path.join(ribo_db_dir, file)
        with open(path, 'r', encoding='utf-8') as ribo_file:
            ribo_elts = parser(ribo_file.read())
        os.remove(path)

        if ribo_elts == "Timeout":
            continue

        out_path = os.path.join(ribo_db_dir, file.rsplit('.', 1)[0] + '.fasta')
        save_to_fasta(ribo_elts, out_path)
        progressbar(n, len(files))

    print("\n   Ribosome database completed\n")
    return n
            
#prepare_ribo_db(r"C:\Subpbiotech_cours\BT5\BIM_BMC\GENOM\project\project_git\GENOM\input\ribosome_db\Brinkman")

