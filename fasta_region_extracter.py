import os
import gzip
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm


class FastaExtract(object):
    def __init__(self, args):
        args = args
        self.mafft_alignment = args.alignment
        self.start_position = args.start
        self.end_position = args.end
        self.cutoff = args.cutoff
        self.ref = args.reference

        # Data
        self.roi_dict = dict

        # Output files
        self.extracted = '.'.join(self.mafft_alignment.split('.')[:-1]) + '_extracted.txt'
        self.stats = os.path.dirname(self.mafft_alignment) + '/ROI_stats.tsv'
        self.output_tsv = os.path.dirname(self.mafft_alignment) + '/ROI_stats.tsv'  # Create path to report

        # Run
        self.run()

    def run(self):
        self.check()  # Check if can find input alignment file

        # Parse alignment file to extract counts for each variant of the region of interest
        print('Parsing alignment file and extraction ROI...')
        FastaExtract.extract(self.mafft_alignment, self.start_position, self.end_position)
        roi_dict = FastaExtract.filter_n(self.extracted)

        # Output the variant frequency table
        print('Filtering ROI and preparing report file...')
        df1 = FastaExtract.filter_cutoff(roi_dict, self.cutoff, self.ref)

        # Replace conserved bases with dots
        roi_ref = FastaExtract.extract_ref(self.ref, self.start_position, self.end_position)
        df2 = FastaExtract.mask_alignment(df1, roi_ref)

        # Print table
        FastaExtract.print_table(df2, self.output_tsv)

    def check(self):
        if '~' in self.mafft_alignment:
            self.mafft_alignment = os.path.expanduser(self.mafft_alignment)
        if not os.path.isfile(self.mafft_alignment):
            raise Exception('Could not locate alignment file. Please use absolute path to file.')

    @staticmethod
    def block_read(my_file, size=1048576):
        while True:
            b = my_file.read(size)
            if not b:
                break
            yield b

    @staticmethod
    def count_lines(my_file):
        total_cnt = 0  # Total entries in fasta file
        with gzip.open(my_file, 'rb', 1024*1024) if my_file.endswith('gz') \
                else open(my_file, 'r', 1024*1024) as in_f:
            if my_file.endswith('gz'):
                total_cnt = sum(bl.count(b'\n') for bl in FastaExtract.block_read(in_f))
            else:
                total_cnt = sum(bl.count('\n') for bl in FastaExtract.block_read(in_f))
        return total_cnt

    @staticmethod
    def extract(input_alignment, start, end):
        out_file = '.'.join(input_alignment.split('.')[:-1]) + '_extracted.txt'

        start_index = start - 1
        end_index = end - 1

        with open(out_file, 'w') as out_f:
            with gzip.open(input_alignment, 'rb', 1024*1024) if input_alignment.endswith('gz') \
                    else open(input_alignment, 'r', 1024*1024) as in_f:
                seq_list = list()
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
                        out_f.write('{}\n'.format(extracted_seq))
                        seq_list = list()  # Empty list that collects the
                    else:  # Sequence line
                        seq_list.append(line)  # Fetch sequences spanning over multiple lines
                # For the last entry in the file
                seq = ''.join(seq_list)
                extracted_seq = seq[start_index:end_index].upper()
                out_f.write('{}\n'.format(extracted_seq))

    @staticmethod
    def filter_n(extracted_roi):
        seq_dict = dict()

        with open(extracted_roi, 'r', 1024*1024) as in_f:
            for line in in_f:
                line = line.rstrip()
                if not line:  # If line is empty
                    continue  # Skip to next line
                if any([x in line for x in ['N', 'n']]):  # Skip extracted sequences with Ns
                    continue

                # Add to dictionary
                if line in seq_dict:
                    seq_dict[line] += 1  # Increase count by 1
                else:
                    seq_dict[line] = 1  # Set count to 1

        return seq_dict

    @staticmethod
    def extract_ref(ref, start, end):
        with gzip.open(ref, 'rb', 1024*1024) if ref.endswith('gz') \
                else open(ref, 'r', 1024*1024) as in_f:
            seq_list = list()
            for line in in_f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    continue
                else:
                    seq_list.append(line)
        return ''.join(seq_list)[start - 1: end - 1]

    @staticmethod
    def filter_cutoff(seq_dict, cutoff, roi_ref):
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

        # Insert reference sequence at first line
        old_idx = df.index
        new_idx = ['Wuhan', old_idx]
        df.loc[-1] = ['', '', roi_ref]  # adding a row

        # with open(output_tsv, 'w') as out_f:
        #    df.to_csv(out_f, sep='\t', header=True, index=True)
        return df

    @staticmethod
    def mask_alignment(df):
        # compare to wild type
        return df

    @staticmethod
    def print_table(df, output_tsv):
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
    parser.add_argument('-r', '--reference', metavar='reference.fasta',
                        required=True,
                        type=str,
                        help='Reference fasta file.'
                             ' Mandatory.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FastaExtract(arguments)
