from Bio import SeqIO
import gzip
import glob
import argparse
import re
import csv
import os
from collections import defaultdict
from tqdm import tqdm
from functools import lru_cache

parser = argparse.ArgumentParser(description='Demultiplex fasta based on the header. '
                                             'This version uses i7 index from paired barcode')
parser.add_argument('-f', '--fastq_folder', default='.', help='Folder containing fastq.gz files for demultiplex')
parser.add_argument('-o', '--output', help='Output folder')
parser.add_argument('-b', '--barcodes', help='File with barcodes that should be demultiplexed and samle names. '
                                             'Other barcodes will be discarded')
parser.add_argument('-d', '--distance', default=0, type=int, help='Distance between observed and target barcode')
args = parser.parse_args()

# Create dir for demuxed reads
if not os.path.exists(args.output):
    os.makedirs(args.output)

# read barcodes and sample names for demultiplexing
barcodes = defaultdict(str)
barcodes_list = list()
with open(args.barcodes, 'r') as b_file:
    barcodes_file = csv.reader(b_file, delimiter='\t')
    for bc in barcodes_file:
        barcodes[bc[0]] = bc[1]
        barcodes_list.append(bc[0])


# distance calculation function one-to-many
@lru_cache(maxsize=100000)
def bc_distance(bc1):
    bc_distance = [sum([1 for x, y in zip(bc1, bcl) if x.lower() != y.lower()]) for bcl in barcodes.keys()]
    return min(bc_distance)


# collect fastq files
fastq_files = glob.glob(args.fastq_folder + '/*.fq.gz')

for fq_file in fastq_files:
    # extract lane number from fastq file name
    match = re.search(r'_(L\d+)_(\d+)', fq_file)
    lane = match.group(1)
    pair = match.group(2)

    with gzip.open(fq_file, 'rt') as fq_file_open:
        print(f'Analyzing {fq_file}')

        n_reads = int(int(os.popen('zcat < ' + fq_file + ' | wc -l').read()) / 4)
        count_demux_reads = 0

        gz_out = dict()
        for barcode in barcodes:
            gz_out[barcode] = gzip.open(args.output + '/' + barcode + '_' + lane + '_' + pair + '.fq.gz', 'at')

        undet_hndl = gzip.open(args.output + '/' + 'Undetermined_' + lane + '_' + pair + '.fq.gz', 'at')

        with tqdm(total=n_reads) as pbar:
            for record in SeqIO.parse(fq_file_open, "fastq"):
                pbar.update(1)
                # Collect i7 read
                read_i7 = record.description.split(':')[9].split('+')[0]
                # Check if barcode exists in the sample set
                if args.distance:
                    distance = bc_distance(read_i7)
                else:
                    if not barcodes.get(read_i7) is None:
                        distance = 0
                    else:
                        distance = 100

                if distance <= args.distance:
                    # count demultiplexed reads
                    count_demux_reads += 1
                    # Write barcode to the separate file
                    SeqIO.write(sequences=[record], handle=gz_out[read_i7], format="fastq")
                else:
                    # Write non-demultiplexed barcodes to separate file
                    SeqIO.write(sequences=[record], handle=undet_hndl, format="fastq")

        for barcode in barcodes:
            gz_out[barcode].close()

        undet_hndl.close()

        p_demux = count_demux_reads / n_reads * 100
        print(f'Demultiplexed {p_demux:.2f}% of reads in {fq_file}')
