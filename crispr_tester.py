import os
import gzip
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
import subprocess
import io
from math import floor


class FastaExtract(object):
    def __init__(self, args):
        args = args
        self.mafft_alignment = args.alignment
        self.query = args.query
        self.cutoff = args.cutoff
        self.ref = args.reference
        self.gap = args.gap
        self.mismatch = args.mismatch

        # Data
        self.roi_dict = dict

        # Output files
        self.extracted = '.'.join(self.mafft_alignment.split('.')[:-1]) + '_extracted.txt'
        self.stats = os.path.dirname(self.mafft_alignment) + '/ROI_stats.tsv'
        self.output_tsv = os.path.dirname(self.mafft_alignment) + '/ROI_stats.tsv'
        self.masked_output_tsv = os.path.dirname(self.mafft_alignment) + '/ROI_masked_stats.tsv'
        self.gggenome_raw_output_tsv = os.path.dirname(self.mafft_alignment) + '/GGGenome_matches.tsv'
        self.gggenome_summary_output_tsv = os.path.dirname(self.mafft_alignment) + '/GGGenome_summary.tsv'

        # Run
        self.run()

    def run(self):
        self.check()  # Check if can find input alignment file

        # Compute the maximum number of mismatch(es) automatically based on query length
        max_diff = int(floor(len(self.query) * self.mismatch))

        ##### Cross-reactivity test #####

        print('Performing cross-reactivity test for {}'.format(self.query))

        # GGGenome
        ggg_df, diff_dict = FastaExtract.loop_gggenome(self.query, max_diff, self.gap)

        # Write GGGenome df to CSV file
        FastaExtract.print_table(ggg_df, self.gggenome_raw_output_tsv)

        # Create GGGenome summary dataframe
        sum_df = FastaExtract.create_gggenome_summary(diff_dict, max_diff)

        # Write GGGenome summary to file from dataframe
        FastaExtract.print_table(sum_df, self.gggenome_summary_output_tsv)

        # Extract start and end positions from the reference Wuhan sequences to perform the inclusivity analysis
        start, end = FastaExtract.get_start_end_positions(ggg_df)


        ##### Inclusivity test #####

        print('\nPerforming inclusivity test')

        roi_ref = FastaExtract.extract_ref(self.ref, start, end)
        print('\tReference ROI: {} ({}..{})'.format(roi_ref, start, end))

        # Parse alignment file to extract counts for each variant of the region of interest
        print('\tParsing MAFFT alignment file and extracting ROI for each entry...')
        # FastaExtract.extract(self.mafft_alignment, self.start_position, self.end_position)  # comment for debug
        roi_dict = FastaExtract.filter_n(self.extracted)

        # Output the variant frequency table
        print('\tFiltering ROI and preparing report file...')
        # Replace conserved bases with dots
        df = FastaExtract.filter_cutoff(roi_dict, self.cutoff, roi_ref)
        FastaExtract.print_table(df, self.output_tsv)  # Print table

        # Mask sequences
        FastaExtract.mask_alignment(df)
        FastaExtract.print_table(df, self.masked_output_tsv)

    def check(self):
        if '~' in self.mafft_alignment:
            self.mafft_alignment = os.path.expanduser(self.mafft_alignment)
        if not os.path.isfile(self.mafft_alignment):
            raise Exception('Could not locate alignment file. Please use absolute path to file.')

        # Check if mismatches/gaps is 25% or less. Otherwise GGGenome will return an error.
        if self.mismatch > 0.25:
            raise Exception('Number of mismatches/gaps should be 25% (0.25) or less.')

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
        extracted_ref = ''.join(seq_list)[start - 1: end - 1]
        return extracted_ref

    @staticmethod
    def filter_cutoff(seq_dict, cutoff, roi_ref):
        # Convert dictionary to Pandas dataframe
        df = pd.DataFrame.from_dict(seq_dict, orient='index')

        # Add column at the end for frequencies
        total = df[0].sum()  # Sum all the counts
        freq_list = [float('{:.2f}'.format(x/total*100)) for x in df[0]]
        df['Frequency (%)'] = freq_list

        # Compute the added % of the lines to remove (below the cutoff)
        df_rem = df[df['Frequency (%)'] < cutoff * 100]
        cnt_rem = df_rem[0].sum()
        freq_rem = float('{:.2f}'.format(cnt_rem/total*100))

        # Remove lines with frequency below cutoff
        df = df[df['Frequency (%)'] >= cutoff * 100]

        # Sort descending based on index values
        df.sort_values(by=0, ascending=False, inplace=True)

        # Add a new index column
        df.reset_index(inplace=True)

        # Change index values
        df.index = np.arange(1, len(df[0]) + 1)

        # Rename "index" column
        df = df.rename(columns={'index': 'Variant'})

        # Rename columns
        df.columns = ['Variant', 'Count', 'Frequency (%)']

        # Add name to index column
        df.index.rename('Group', inplace=True)

        # Move columns
        df = df[['Count', 'Frequency (%)', 'Variant']]

        # Insert reference sequence at first line
        new_row = pd.DataFrame({'Count': '', 'Frequency (%)': '', 'Variant': roi_ref}, index=['Wuhan-Hu-1'])
        df = pd.concat([new_row, df])

        # Add a row at the end named "Removed" with the added % of the removed lines so the added freq.=100%
        rem_row = pd.DataFrame({'Count': cnt_rem, 'Frequency (%)': freq_rem, 'Variant': ''},
                               index=['Variants with frequency lower than {}%'.format(cutoff * 100)])
        df = pd.concat([df, rem_row])

        return df

    @staticmethod
    def mask_alignment(df):
        # Get ref sequence
        ref = df.loc['Wuhan-Hu-1', 'Variant']
        # compare to wild type
        for i, row in df.iterrows():
            if i == 'Wuhan-Hu-1':  # Don't change the reference sequence
                continue
            seq = df.loc[i, 'Variant']
            seq_list = list()
            for j, b in enumerate(seq):
                if b == ref[j]:
                    seq_list.append('.')
                else:
                    seq_list.append(b)
            masked_seq = ''.join(seq_list)

            # Replace variant sequence in dataframe with masked one
            df.loc[i, 'Variant'] = masked_seq

        # print(df)
        return df

    @staticmethod
    def print_table(df, output_tsv):
        with open(output_tsv, 'w') as out_f:
            df.to_csv(out_f, sep='\t', header=True, index=True)

    @staticmethod
    def run_mafft(input_fasta, reference_fasta, mafft_cutoff):
        output_alignment = '.'.join(input_fasta.split('.')[:-1]) + '_mafft' + mafft_cutoff + '.fasta'
        cmd = ['mafft', '--auto', '--keeplength',
               '--maxambiguous', str(mafft_cutoff),
               '--addfragments', input_fasta, reference_fasta]
        with open(output_alignment, 'w') as f:
            subprocess.Popen(cmd, stdout=f)  # write standard output (alignment) to file

    @staticmethod
    def run_gggenome_online(seq, db, mismatch, gap=False):
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
        # http://gggenome.dbcls.jp/COVID19-primercheck-EUA-20200501/5/TTTGCCCCCAGCGCTTCAGCGTT
        # http://gggenome.dbcls.jp/hg38/4/nogap/TTTGCCCCCAGCGCTTCAGCGTT
        url = 'https://GGGenome.dbcls.jp/{}/{}/nogap/{}.csv.download'.format(db, mismatch, seq)
        if gap:
            url = 'https://GGGenome.dbcls.jp/{}/{}/{}.csv.download'.format(db, mismatch, seq)
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            r.raise_for_status()
            # raise Exception('Problem with GGGenome URL request: {}'.format(r.status_code))
        else:
            # Parse results into Pandas dataframe using "fake" file handle with SingIO
            try:
                df = pd.read_csv(io.StringIO(r.content.decode()), sep=',', skiprows=6)
                return df
            except Exception as e:
                print(type(e).__name__, e)

    @staticmethod
    def loop_gggenome(roi_ref, max_diff, gap):
        # Loop dataframe and check sequence with GGGenome
        # Results are returned in a new dataframe

        # Master GGGenome dataframe
        ggg_df = pd.DataFrame(columns=['# name', 'strand', 'start', 'end', 'snippet', 'snippet_pos', 'snippet_end',
                                       'query', 'sbjct', 'align', 'edit', 'match', 'mis', 'del', 'ins'])

        # Reformat GGGenome output table
        # Loop through GGGenome results for each query with increasing differences allowed
        # TODO: check for hits in human genome too.
        db_list = ['COVID19-primercheck-EUA-20200501', 'hg38']
        for db in db_list:
            print('\tTesting query online with GGGenome using "{}" database'.format(db))
            for diff in range(0, max_diff+1):
                df1 = FastaExtract.run_gggenome_online(roi_ref, db, diff, gap)
                # Check if results returned
                if df1.empty or df1['# name'][0] == '### No items found. ###':
                    continue
                else:
                    ggg_df = pd.concat([ggg_df, df1])

                # # Check if df1 is not empty
                # if df1.empty:
                #     raise Exception('Could not find any match in GGGenome "{}" database.'.format(db))
                #
                # # Concatenate with master GGGenome dataframe
                # ggg_df = pd.concat([ggg_df, df1])

        # Remove duplicated entries
        ggg_df = ggg_df.drop_duplicates()

        # Reset index. Some rows have to same index.
        ggg_df.reset_index(drop=True, inplace=True)

        # # Remove matches to reference
        # to_drop_list = ggg_df.index[ggg_df['sbjct'] == (roi_ref or FastaExtract.reverse_complement(roi_ref))].tolist()
        # ggg_df.drop(to_drop_list, inplace=True)
        # ggg_df.reset_index(drop=True, inplace=True)  # Reset pandas index

        # Add results to summary dataframe
        name_list = ggg_df['# name'].to_list()  # convert name column to list
        org_list = [FastaExtract.extract_org(x) for x in name_list]

        # remove duplicates from list
        org_list = list(dict.fromkeys(org_list))
        diff_dict = dict()

        for org in org_list:
            if org not in diff_dict.keys():
                diff_dict[org] = dict()
                # Generate headers for differences (mismatches and gaps)
                for diff in range(0, max_diff+1):
                    diff_dict[org][diff] = 0

            # Fetch the difference values
            matching_row_list = ggg_df.index[ggg_df['# name'].str.contains(org)].tolist()
            for j in matching_row_list:
                mismatches = len(roi_ref) - ggg_df.iloc[j]['match']
                diff_dict[org][mismatches] += 1  # Add 1 to the count

        return ggg_df, diff_dict

    @staticmethod
    def create_gggenome_summary(diff_dict, max_diff):
        # Generate headers for differences (mismatches and gaps)
        header_list = list()
        for diff in range(0, max_diff+1):
            header_list.append(str(diff))

        # Convert dictionary to pandas dataframe
        df = pd.DataFrame.from_dict(diff_dict, orient='index')

        # Rename index column
        df = df.rename_axis("Organism")
        # df.index.rename(['Organism'], inplace=True)

        # Prepare name for columns
        header_list = ['m' + str(x) for x in header_list]

        # Rename columns
        df.columns = header_list

        # Add a column at the end with the total mismatches per organism
        df['Total'] = df.sum(axis=1)

        # Sort rows by 1) total number of mismatches and 2) alphabetically
        df.sort_values(by=header_list, ascending=[False] * len(header_list), inplace=True)
        # df.sort_index(ascending=True, inplace=True)  # Sort alphabetically

        # Add a row at the end with total for each mismatch column
        df.loc['Total'] = df.sum(numeric_only=True, axis=0)

        return df

    @staticmethod
    def reverse_complement(seq):
        rc_tab = str.maketrans('ATCG', 'TAGC')
        return seq.translate(rc_tab)[::-1]

    @staticmethod
    def extract_org(header):
        # regex = re.compile(r"\s*ctg\s*", flags=re.IGNORECASE)
        org = header.split(' ')
        if len(org) > 1:
            org = header.split(' ', 1)[1]  # Ditch the accession number
            org = org.split(',')[0]
            org = org.split('chromosome')[0]
            org = org.split('NODE')[0]
            org_list = org.split()  # split string into list
            for i, word in enumerate(org_list):
                if any([substring in word for substring in ['ctg', 'scaf', 'scf']]):
                    del org_list[i]
            org = ' '.join(org_list)
            org = org.split('contig')[0]
            return org
        else:  # hit in human genome (eg. chr2)
            return header

    @staticmethod
    def extract_acc(header):
        return header.split(' ')[0]

    @staticmethod
    def get_start_end_positions(df):
        return df.loc[0, 'start'], df.loc[0, 'end']


if __name__ == '__main__':
    parser = ArgumentParser(description='Extract region from fasta file.')
    parser.add_argument('-a', '--alignment', metavar='mafft_alignment.fasta',
                        required=True,
                        type=str,
                        help='Mafft alignment file.'
                             ' Mandatory.')
    parser.add_argument('-q', '--query', metavar='TTTNCCCCCAGCGCTTCAGCGTTC',
                        required=True,
                        type=str,
                        help='Sequence to test (PAM+crRNA).'
                             ' Mandatory.')
    parser.add_argument('-c', '--cutoff', metavar='0.001',
                        type=float, default=0.001,
                        required=False,
                        help='Cutoff frequency to keep a variant. Must be between 0 and 1.'
                             ' Default is 0.001.')
    parser.add_argument('-r', '--reference', metavar='reference.fasta',
                        required=True,
                        type=str,
                        help='Reference fasta file.'
                             ' Mandatory.')
    parser.add_argument('-g', '--gap',
                        action='store_true',
                        help='Allows gaps in GGGenome search. Default is False.')
    parser.add_argument('-m', '--mismatch', metavar='0.20',
                        type=float,
                        required=False,
                        default=0.20,
                        help='Percentage of mismatch allowed.'
                             ' Minimum is 0 and maximum is 0.25 (GGGenome requirement).'
                             ' Includes gaps if "-g" is used.'
                             ' Default is 0.20.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    FastaExtract(arguments)
