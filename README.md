# CRISPR_Assay_Tester
Test CRISPR assay specificity and sensitivity on small genomes

## Description
This is a python3 script.
Use mafft(https://mafft.cbrc.jp/alignment/software/), gggenome (https://gggenome.dbcls.jp/)

## Installation
- Make sure you have `conda` installed (https://docs.conda.io/en/latest/miniconda.html).
- Create a conda virtual environment and install required packages:
```
conda create -n crispr_tester git panadas mafft requests tqdm
```
- Clone repository and test tool:
```
git clone https://github.com/duceppemo/CRISPR_Assay_Tester
cd CRISPR_Assay_Tester
python crispr_tester.py -h
```

## Usage
```
usage: crispr_tester.py [-h] -a mafft_alignment.fasta -s 29192 -e 29215 -c 0.01 -r reference.fasta

Extract region from fasta file.

optional arguments:
  -h, --help            show this help message and exit
  -a mafft_alignment.fasta, --alignment mafft_alignment.fasta
                        Mafft alignment file. Mandatory.
  -s 29192, --start 29192
                        Start position for extraction. Mandatory.
  -e 29215, --end 29215
                        Length of the sequence to extract. Mandatory.
  -c 0.01, --cutoff 0.01
                        Cutoff frequency to keep a variant. Must be between 0 and 1. Default is 0.01. Mandatory.
  -r reference.fasta, --reference reference.fasta
                        Reference fasta file. Mandatory.
```

## Outputs

Please set the font to "Courier new" ((monospaced font) to view the tables for proper sequenced alignment display.
