from ete3 import NCBITaxa
from taxoniq import Taxon

def load_taxname(path : str) -> list:
    with open(path, 'r') as taxfile:
        taxlist = [line.strip('\n') for line in taxfile.readlines()]

    return taxlist

def get_id(path : str):
    taxlist = load_taxname(path)
    return NCBITaxa().get_name_translator(taxlist)

def get_accession(taxid : int):
    return Taxon(taxid).refseq_genome_accessions[0].accession_id

