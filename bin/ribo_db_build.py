import re, os
from dataclasses import dataclass
from typing import List, Dict, Tuple, Generator

@dataclass
class GenBankElement:
    organism_name : str
    accession_id : str
    ribosome_beacon : Dict[str, Tuple[int, int]]
    ribosome_sequences : Dict[str, str]
    
def load_files(ribo_db_dir : str) -> Generator:
    """
    Load one after the other the GenBank general files to parse.
    """
    for file in os.listdir(ribo_db_dir):
        if file.endswith('.gb'):
            filepath = os.path.join(ribo_db_dir, file)
            with open(filepath, 'r', encoding='utf-8') as ribo_file:
                yield ribo_file.read()
        os.remove(filepath)

def postprocess_regex(pattern : List[Tuple[str]]) -> List[List[str]]:
    liste = []
    for elt in pattern:
        new_tup = [i for i in elt if i != '']
        liste.append(new_tup)
    return liste

def parser(ribo_db_dir : str):
    """
    Parses the .gb files downloaded with the other genomes and
    extract the ribosomic DNA sequences along with other informations.
    Uses RegEx for parsing : sensible to GenBank Summary format.
    """
    for ribo_file_content in load_files(ribo_db_dir):
        genome_seq = ''.join(re.findall('[atcg]+', ribo_file_content.split('ORIGIN')[1]))
        accession_id = re.findall(r'VERSION     (.*?(?=\n))', ribo_file_content)[0]
        organism_name = re.findall(r'ORGANISM..(.+?)(?=\n)', ribo_file_content)[0]
        pattern = re.findall(r'(?s)rRNA\s+(?:complement\(join\(([\d\.\.,\s]+)\)\)|complement\((\d+\.\.\d+)\)|(\d+\.\.\d+)).*?product=\"(.*?)(?= ribosomal RNA)',
                             ribo_file_content)
        pattern = postprocess_regex(pattern)
        ribo_dico = {}
        for elt in pattern:
            beacons, ribo_id = elt
            beacons_list = re.findall('([0-9]+..[0-9]+)', beacons)

            ribo_dico[ribo_id] = [(int(beacons.split('..')[0]), int(beacons.split('..')[1])) for beacons in beacons_list]
        
        ribo_seq_dico = {ribo_id : [genome_seq[beacon[0]:beacon[1]] for beacon in beacons] for ribo_id, beacons in ribo_dico.items()}      
        
        yield GenBankElement(
                ribosome_sequences=ribo_seq_dico,
                organism_name=organism_name,
                ribosome_beacon=ribo_dico,
                accession_id=accession_id                
            )

def prepare_ribo_db(ribo_db_dir : str):
    for ribo_elts in parser(ribo_db_dir):
        name_file = f'ribosomes_{ribo_elts.organism_name}.fasta'
        out_path = os.path.join(ribo_db_dir, name_file)
        if name_file in os.listdir(ribo_db_dir):
            os.remove(out_path)
        
        for ribo_id, beacons in ribo_elts.ribosome_beacon.items():
            for index, seq in enumerate(ribo_elts.ribosome_sequences[ribo_id]):
                with open(out_path, 'a+') as ribofile:
                    ribofile.write(
                        f">{ribo_id}_{index}_{ribo_elts.organism_name}|{beacons[index][0]}-{beacons[index][1]}\n{seq.upper()}\n"
                    )
    return 0
            



