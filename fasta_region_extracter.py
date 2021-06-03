import os
import gzip
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
import mmap


class FastaExtract(object):
    def __init__(self, args):
        args = args
        self.mafft_alignment = args.alignment
        self.start_position = args.start
        self.end_position = args.end
        self.cutoff = args.cutoff

        # Data
        self.roi_dict = dict

        # Run
        self.run()

    def run(self):
        self.check()  # Check if can find input alignment file
        output_tsv = os.path.dirname(self.mafft_alignment) + '/ROI_stats.tsv'  # Create path to report

        # Parse alignment file to extract counts for each variant of the region of interest
        print('Parsing alignment file and extraction ROI...')
        roi_dict = self.extract(self.mafft_alignment, self.start_position, self.end_position)

        # Output the variant frequency table
        print('Filtering ROI and preparing report file...')
        self.report(roi_dict, output_tsv, self.cutoff)

    def check(self):
        if '~' in self.mafft_alignment:
            self.mafft_alignment = os.path.expanduser(self.mafft_alignment)
        if not os.path.isfile(self.mafft_alignment):
            raise Exception('Could not locate alignment file. Please use absolute path to file.')

    @staticmethod
    def block_read(my_file, size=65536):
        while True:
            b = my_file.read(size)
            if not b:
                break
            yield b

    @staticmethod
    def count_fasta_entries(fasta_file):
        total_cnt = 0  # Total entries in fasta file
        with gzip.open(fasta_file, 'rb') if fasta_file.endswith('gz') else open(fasta_file, 'r') as in_f:
            if fasta_file.endswith('gz'):
                total_cnt = sum(bl.count(b'>') for bl in FastaExtract.block_read(in_f))
            else:
                total_cnt = sum(bl.count('>') for bl in FastaExtract.block_read(in_f))
        return total_cnt

    @staticmethod
    def count_lines(my_file):
        total_cnt = 0  # Total entries in fasta file
        with gzip.open(my_file, 'rb') if my_file.endswith('gz') else open(my_file, 'r') as in_f:
            if my_file.endswith('gz'):
                total_cnt = sum(bl.count(b'\n') for bl in FastaExtract.block_read(in_f))
            else:
                total_cnt = sum(bl.count('\n') for bl in FastaExtract.block_read(in_f))
        return total_cnt

    @staticmethod
    def get_num_lines(file_path):
        fp = open(file_path, "r+")
        buf = mmap.mmap(fp.fileno(), 0)
        lines = 0
        while buf.readline():
            lines += 1
        return lines

    @staticmethod
    def extract(input_alignment, start, end):
        seq_dict = dict()  #

        start_index = start - 1
        end_index = end - 1

        print('Conting entries...')
        #total_cnt = FastaExtract.count_fasta_entries(input_alignment)

        print('Parsing file...')
        with gzip.open(input_alignment, 'rb') if input_alignment.endswith('gz') else open(input_alignment, 'r') as in_f:
            seq_list = list()
            #for line in in_f:  # Step input file line by line
            #for line in tqdm(in_f, total=FastaExtract.get_num_lines(input_alignment)):
            for line in tqdm(in_f, total=FastaExtract.count_lines(input_alignment)):
                line = line.rstrip()
                if input_alignment.endswith('gz'):
                    line = line.decode()  # Convert from binary to string
                if not line:  # If line is empty
                    continue  # Skip to next line

                if line.startswith('>') and not seq_list:  # First line
                    continue  # Skip to next line
                elif line.startswith('>') and seq_list:  # A new sequence starts
                    seq = ''.join(seq_list)  # Combine all lines if sequence spanned over multiple lines
                    extracted_seq = seq[start_index:end_index].upper()  # Extract region of interest by slicing
                    if 'N' in extracted_seq:  # Skip extracted sequences with Ns
                        continue

                    # Add to dictionary
                    try:  # Try to add it
                        seq_dict[extracted_seq] += 1  # Increase count by 1
                    except KeyError:  # Create key if it doesn't exist
                        seq_dict[extracted_seq] = 1  # Set count to 1

                    seq_list = list()  # Empty list that collects the
                else:  # Sequence line
                    seq_list.append(line)  # Fetch sequences spanning over multiple lines
            # For the last entry in the file
            seq = ''.join(seq_list)
            extracted_seq = seq[start_index:end_index].upper()
            if 'N' in extracted_seq:  # Skip extracted sequences with Ns
                return  # ???
            try:
                seq_dict[extracted_seq] += 1
            except KeyError:
                seq_dict[extracted_seq] = 1
        return seq_dict

    @staticmethod
    def report(seq_dict, output_tsv, cutoff):
        # Convert dictionary to Pandas dataframe
        df = pd.DataFrame.from_dict(seq_dict, orient='index')

        # Add column at the end for frequencies
        total = df[0].sum()  # Sum all the counts
        freq_list = [float('{:.2f}'.format(x/total*100)) for x in df[0]]
        df['Frequency'] = freq_list

        # Remove lines with frequency below cutoff
        df = df[df['Frequency'] > cutoff * 100]

        # Sort descending based on index values
        df.sort_values(by=0, ascending=False, inplace=True)

        # Add a new index column
        df.reset_index(inplace=True)

        # Change index values
        df.index = np.arange(1, len(df[0]) + 1)

        # Add name to index column
        df.index.rename('Group', inplace=True)

        # Rename "index" column
        df = df.rename(columns={'index': 'Variant'})

        # Rename columns
        df.columns = ['Variant', 'Count', 'Frequency']

        # Move columns
        df = df[['Count', 'Frequency', 'Variant']]
        print(df)

        with open(output_tsv, 'w') as out_f:
            df.to_csv(out_f, sep='\t', header=True, index=True)

    @staticmethod
    def run_gggenome_online(seq, output_file):
        """
        http://gggenome.dbcls.jp/help.html

        http[s]://GGGenome.dbcls.jp/db/k/[strand]/[nogap]/sequence[.format][.download]
            db: hg19, mm10, dm3, ce10, TAIR10, pombe, refseq, etc. (default: hg19)
            Full list of databases is available below.
            k: Maximum number of mismatches/gaps. (default: 0)
            strand: '+' ('plus') or '-' ('minus') to search specified strand only. (optional)
            nogap: No gaps. Allow only k mismatches. (optional)
            sequence: Nucleotide sequence, case insensitive.
            format: html, txt, csv, bed, gff, json. (default: html)
            download: Download result as a file. (optional)
        """
        url = 'http[s]://GGGenome.dbcls.jp/SARS-CoV-2/0/+/nogap/{}.csv.download'.format(seq)
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            r.raise_for_status()
            # raise Exception('Problem with GGGenome URL request: {}'.format(r.status_code))
        else:
            with open(output_file, 'w') as f:
                f.write(r.content.decode())


if __name__ == '__main__':
    parser = ArgumentParser(description='Extract region from fasta file.')
    parser.add_argument('-a', '--alignment', metavar='mafft_alignment.fasta',
                        required=True,
                        type=str,
                        help='Mafft alignment file.'
                             ' Mandatory.')
    parser.add_argument('-s', '--start', metavar='23',
                        required=True,
                        type=int,
                        help='Start position for extraction.'
                             ' Mandatory.')
    parser.add_argument('-e', '--end', metavar='100',
                        required=True,
                        type=int,
                        help='Length of the sequence to extract.'
                             ' Mandatory.')
    parser.add_argument('-c', '--cutoff', metavar=0.01,
                        type=float, default=0.01,
                        required=True,
                        help='Cutoff frequency to keep a variant. Must be between 0 and 1. Default is 0.01 (1%).'
                             ' Mandatory.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FastaExtract(arguments)
