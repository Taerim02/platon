import argparse
import csv
import os

import pyfastx

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("contig_file", help="Input FASTA file")
parser.add_argument("ref_output", help="ref ouput file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()

contigs= {}
tsv_header = ['contig', 'contig_start', 'contig_end', 'plasmid_start', 'plasmid_end', 'plasmid_id', 'plasmid_length', 'coverage', 'identity']

for record in pyfastx.Fasta(args.contig_file):
    length = len(record.seq)
    contig = {
        'id': record.name,
        'length': length}
    contigs[record.name] = contig

if os.path.getsize(args.ref_output) > 0:
    with open(args.ref_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            contig_id = cols[0]
            if contigs.get(contig_id, False):
                hit = {
                    'contig_start': int(cols[2]),
                    'contig_end': int(cols[3]),
                    'plasmid_start': int(cols[4]),
                    'plasmid_end': int(cols[5]),
                    'plasmid_id': cols[1],
                    'plasmid_length': int(cols[6]),
                    'coverage': float(cols[7]) / contigs[contig_id]['length'],
                    'identity': float(cols[8]) / float(cols[7])
                }
                if(hit['coverage'] >= 0.8):
                    with open(args.output, "a", newline='') as fh:
                        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                        is_empty = fh.tell() == 0
                        if is_empty:
                            writer.writeheader()
                        tsv_row = {"contig":contig_id, 'contig_start': hit['contig_start'], 'contig_end':hit['contig_end'], 'plasmid_start': hit['plasmid_start'], 
                            'plasmid_end': hit['plasmid_end'], 'plasmid_id':hit['plasmid_id'], 'plasmid_length':hit['plasmid_length'],
                                    'identity': hit['identity'], 'coverage':hit['coverage']}
                        writer.writerow(tsv_row)
                
                    
 
