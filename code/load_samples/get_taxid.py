from ete3 import NCBITaxa
from taxoniq import Taxon
import sqlite3

import re, os

from typing import Dict, List

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database("taxdump.tar.gz")

def file_to_list(path : str) -> List[str]:
    """ Load the content of a file into a list of strings. """
    with open(path, 'r') as file:
        result = [line.strip('\n') for line in file.readlines()]
    return result

def is_seqid(seqid : str) -> bool:
    """ Check if the input is formatted as a sequence accesion id. """
    accession_patterns = {
        'GenBank': r'^[A-Z]{1,2}\d{5,6}$',
        'RefSeq': r'^[NX][CGMRTW]_\d{6,}(\.\d+)?$',
        'WGS': r'^[A-Z]{4}\d{8,}$',
        'SRA': r'^[SED]RR\d{6,}$',
        'GEO': r'^GSM\d{6,}$'
    }
    for pattern in accession_patterns.values():
        if re.match(pattern, seqid):
            return True
    return False


def parse_file(path : str) -> Dict[str, List[str]]:
    """ Parse the file to get the corresponding sequence accession ids. """
    id_list = file_to_list(path)

    # file contains sequence accession ids
    if is_seqid(id_list[0]):
        return {id: id for id in id_list}, []

    # file contains taxon names
    else:
        taxids = ncbi.get_name_translator(id_list)
        seqids = {}
        failed = []
        for taxname in id_list:
            # get the taxid
            try:
                taxid = taxids[taxname][0]
            except KeyError as e:
                try : taxid = ncbi.get_fuzzy_name_translation(taxname)[0]
                except sqlite3.OperationalError as e:
                    # by default, fuzzy search is not compiled when downloading the ete3 package
                    error_message = str(e)
                    path = error_message.split(":", 1)[0]
                    path = os.path.dirname(path)
                    print(f"To enable fuzzy taxname search, you need to go to the following directory and do the command 'make':\n{path}")
            # get the sequence accession id
            try: seqids[taxname.replace('/', '-')] = Taxon(taxid).refseq_genome_accessions[0].accession_id
            except: failed.append(taxname)
        return seqids, failed