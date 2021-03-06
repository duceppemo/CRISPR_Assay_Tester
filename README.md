# CRISPR_Assay_Tester
CRISPR cross-reactivity and inclusivity *in silico* test.

## Description
This is a python3 (v3.9 used) script meant to be run in a Unix environment. I will not run on a Windows system, unless you run the Ubuntu terminal in Win10 (https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview).

A working internet access is required to run GGGenome (https://gggenome.dbcls.jp/).

In its current implementation, user needs to supply the MAFFT alignment because it can require hundreds of GB of memory to compute. A High Performace Computer (HPC) is likely needed to perform this step successfully.

Use MAFFT (v7.480 used) (https://mafft.cbrc.jp/alignment/software/) with the following command to generate the required alignment file:

```
mafft \
    --auto \
    --keepLength \
    --threads 24 \
    --maxambiguous 0.01 \
    --addfragments fasta_folder/ reference.fasta \
    > mafft_alignment_0.01.fasta
```

A pre-computed alignment file for all SARS-CoV-2 complete sequences from NCBI (*time stamp*) will be available soon for download.

The SARS-CoV-2 Wuhan-Hu-1 reference sequence can be downloaded by clicking this link: [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta).

## Installation
- Make sure you have `conda` installed (https://docs.conda.io/en/latest/miniconda.html).
- Create a virtual environment to install the required packages:
```
conda create -n crispr_tester git pandas mafft requests tqdm
conda activate crispr_tester
```
- Clone repository and test tool:
```
git clone https://github.com/duceppemo/CRISPR_Assay_Tester
cd CRISPR_Assay_Tester
python crispr_tester.py -h
```

## Usage

```
usage: crispr_tester.py [-h] -a mafft_alignment.fasta [-s TTTNCCCCCAGCGCTTCAGCGTTC] [-l sequence_list.tsv] [-c 0.001] -r reference.fasta -o /output/folder [-g] [-m 0.20]

Extract region from fasta file.

optional arguments:
  -h, --help            show this help message and exit
  -a mafft_alignment.fasta, --alignment mafft_alignment.fasta
                        Mafft alignment file. Mandatory.
  -s TTTNCCCCCAGCGCTTCAGCGTTC, --sequence TTTNCCCCCAGCGCTTCAGCGTTC
                        Sequence to test (PAM+crRNA). Must use "-s" or "-l".
  -l sequence_list.tsv, --list sequence_list.tsv
                        A 2-column tab-separated file: Description<tab>Sequence. Must use "-s" or "-l".
  -c 0.001, --cutoff 0.001
                        Cutoff frequency to keep a variant. Must be between 0 and 1. Default is 0.001.
  -r reference.fasta, --reference reference.fasta
                        Reference fasta file. Mandatory.
  -o /output/folder, --output /output/folder
                        Output folder path.
  -g, --gap             Allows gaps in GGGenome search. Default is False.
  -m 0.20, --mismatch 0.20
                        Percentage of mismatch allowed. Minimum is 0 and maximum is 0.25 (GGGenome requirement). Includes gaps if "-g" is used. Default is 0.20.
```
## Input sequence(s)
If you prefer to run the analysis on a single crRNA, you can use the option `-s`. If you want to test multiple crRNAs, you can use the `-l`option. this option will run each sequence automatically one after the other. A parallel implementation might be implemented in the future.

if the `-l` option is selected, the input file for must contain one entry per line. Each line contains two tab-separated columns whenre the 1st column is a description (try to keep it short as it goes intor the file names) and the 2nd column is the sequence. Here an example file (the `<tab>` are just to indicate where the tab character shoul be inserted):
```
Broughton-N<tab>TTTNCCCCCAGCGCUUCAGCGUUC
Broughton-E<tab>GUGGUAUUCUUGCUAGUUAC
Liu-N<tab>CTGCTGCTTGACAGATTGAACCAG
```
## Outputs
Cross-reactivity output files:
- GGGenome_matches.tsv -> GGGenome "raw" results with full details.
- GGGenome_summary.tsv -> GGGenome results summarized by organism and sorted by best match.

Inclusivity output files:
- mafft_alignment_0.01_extracted.txt -> ROI extracted from all the entries in the MAFFT alignment file.
- ROI_stats.tsv -> Frequency of variants, full sequences
- ROI_masked_stats.tsv -> Frequency of variants, masked sequences

Please change the font to "Courier new" (monospaced font) for a proper alignment of the `ROI_masked_stats.tsv`.
