# GENOM
This project investigates the use of kmers to detect horizontal gene transfers.


## Utilisation:

### 1 Setup the working environment

Clone the repository and move in it :

```bash
git clone https://github.com/flopau76/GENOM/
cd GENOM
```

Setup the conda environment:

```bash
conda env create -f genobtain.yaml
conda activate genobtain
```

Add it to your python path:

```bash
cd code
export PYTHONPATH=$(pwd):$PYTHONPATH
```

### 2. Download the required data

The module `load_samples` allows you to download sequences from a list of taxon names or a list of NCBI sequence accesion IDs. For this, run the following command:

```bash
python3 -m load_samples test_samples.txt
```
The argument can either be an absolute path, or the name of a file located in `input/samples_list`. The command downloads the corresponding sequences into `input/sequence_db`. It further creates a report describing the genomes who failed to download. For downstream analysis, the module also downloads the ribosomal DNA of the taxons into `input/ribosomes_db/samples_file/`

A few files are allready provided as an example.

### 3. Generating data with simulated horizontal transfers

For testing purposes, it can be usefull to generate databases with known horizontal gene transfers (HGT). This can be done by:

```bash
python3 -m HGTransfer test_samples
```

By default, this command creates a new folder  `input/sequence_db/generator_db`, which contains the modified sequences and a file `HGT_report.txt` indicating which transfers have been made.

Further arguments allow to specify the probability of horizontal transfers and to enable or not other mutational events.

### 4. Detect potential HGTs with sliding windows

The module `compute_signatures` aims to detect HGT. For this it performs a sliding window over the genomes and compares the signature of the window to the genome wide signature. Windows with an outsanding values are considered as potential HGT

```bash
python3 -m compute_signatures generator_db -k 6 -w 5000 -s 500 -m 1
```

The first, mandatory argument corresponds to the relative path of the sequence folder, starting from `input/sequence_db`.  
`-k`, `-w` and `-s` are optional parameters of the sliding window, respectively the size of the kmers, the size of a window and the step between two windows.  
`-m` indicates wich metrics to use to compare the windows to the genome. See help for which are currently available.  
This will create a pdf with the resuling graph and a json containing the potential HGT, both located in `output/transfer_summary`.

### 5 Analyse the identified HGT and backtrack their origin
Once we have identified horizontal gene transfer, we might be interested in their origin genome. This can be done via the module `analysis`

```bash
python3 -m analysis generator_db_Chi-squared_distance.json generator_db
```

This command takes as an input the previously generated json and the sequence database use to generate it. For each hit from the json, it will look for a potential donnor in among the possible genomes from the database.

By default, the results will be created in `output/analysis`. They contain a text file describing precisely the position of each HGT in donnor and receiver genome, along with a graph summarizing the transfers between the species.
