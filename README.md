# GENOM
using kmers to investigate horizontal transfer


## Utilisation:

### 1 Setup the working environment

Clone the repository and move in it :

```bash
git clone https://github.com/flopau76/GENOM/
cd GENOM
```

Setup the conda environment:

```
conda env create -n genobtain -f genobtain.yaml
conda activate genobtain
```

### 2. Download the required data

The module `load_samples` allows you to download sequences from a list of taxon names with the following command:

```bash
python -m load_samples <samples_file>
```
`<samples_file>` can either be an absolute path, or the name of a file located in `input/samples_list`.
The sequences are donwloaded into `input/sequence_db/samples_file/`. Each taxa has a separate sub-directory with its fasta file. A report is also created describing the genomes who failed to download.
The module also downloads the ribosomal DNA of the taxons into `input/ribosomes_db/samples_file/`

### 3. Generating data with simulated horizontal transfers

For testing purposes, it can be usefull to generate data with known horizontal gen transfers (HGT). This can be done by:

```bash
python -m generate_HGT <database>
```

TODO: complete this. Add info about ground truth.

### 4. Computing kmer-signatures and comparing them

TODO