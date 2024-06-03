import argparse
import csv
import os

"""Search for oriT sequence hits in the metagenomic module of platon."""

parser = argparse.ArgumentParser(description="Process FASTA file for OriT detection")
parser.add_argument("orit_output", help="Input FASTA file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")

args = parser.parse_args()

# set the variable needed to store data
tsv_header = ['contig', 'contig_start', 'contig_end', 'orit_start', 'orit_end', 'orit_id', 'orit_length', 'coverage', 'identity']

# extract information and store into tsv file
if os.path.getsize(args.orit_output) > 0:
    with open(args.orit_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            contig_id = cols[0]
            hit = {
                'contig_start': int(cols[2]),
                'contig_end': int(cols[3]),
                'orit_start': int(cols[4]),
                'orit_end': int(cols[5]),
                'orit_id': cols[1],
                'orit_length': int(cols[6]),
                'coverage': float(cols[7]) / int(cols[6]),
                'identity': float(cols[8]) / float(cols[7])
            }
            if(hit['coverage'] >= 0.9 and hit['identity'] >= 0.9):
                 with open(args.output, "a", newline='') as fh:
                      writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                      is_empty = fh.tell() == 0
                      if is_empty:
                          writer.writeheader()
                      tsv_row = {"contig":contig_id, 'contig_start': hit['contig_start'], 'contig_end':hit['contig_end'], 'orit_start': hit['orit_start'], 
                          'orit_end': hit['orit_end'], 'orit_id':hit['orit_id'], 'orit_length':hit['orit_length'],
                                  'identity': hit['identity'], 'coverage':hit['coverage']}
                      writer.writerow(tsv_row)
              
