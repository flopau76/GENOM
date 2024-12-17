from ete3 import NCBITaxa
from taxoniq import Taxon

from typing import Dict, List

def parse_file(path : str) -> List[str]:
    """ Load the file containing species/strain name. """
    with open(path, 'r') as file:
        lines = [line.strip('\n') for line in file.readlines()]
    return lines

def get_id(path : str) -> Dict[str,str]:
    """ Get the taxid from a species/strain name. """
    taxlist = parse_file(path)
    return NCBITaxa().get_name_translator(taxlist)

def get_accession(taxid : int) -> str:
    """ Return the genome accession ID for a given taxid. """
    return Taxon(taxid).refseq_genome_accessions[0].accession_id

