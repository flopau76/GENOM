# GENOM
using kmers to investigate horizontal transfer


## Utilisation:
### 1. Download the required data

#### Setup the working environment

Clone the repository and move in it :

```
git clone https://github.com/flopau76/GENOM/
cd GENOM
```

Setup the conda environment:

```
conda env create -n genobtain -f genobtain.yaml
conda activate genobtain
```

You can now download the wanted genomes by writing the taxa name in a `.txt` file and entering the following command:

```
#For Linux
./bin/get_data.py ~/PATH/TO/<your_taxa_file>.txt ~/PATH/TO/OUT_DIR
```

In the output directory is created for each taxa a directory containing the fasta file. A report is also outputed describing the number of genomes successfully downloaded and those who failed, indicating which failed due to unavailability or wrong name inputed. An example of input file is given in `./sample/taxa_list.txt`.

### 2. Run the python module
To run the programm on the `./data` folder, directly launch the main module with the following command

```bash
    # Déclenche l'exécution du fichier __main__.py dans le module TP
    python3 -m TP
```
