import os
import pyfastx
import argparse
from pathlib import Path
import subprocess as sp
import sys

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("contig_file", help="Input FASTA file")
parser.add_argument("ref_output", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()
contigs= []
for record in pyfastx.Fasta(args.contig_file):
    length = len(record.seq)
    contig = {
        'id': record.name,
        'length': length}
    contigs.append(contig)
log_file = os.path.join(f'{args.name}.log')
length = 0
if os.path.getsize(args.ref_output) > 0:
    with open(args.ref_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            contig_id = cols[0]
            for contig in contigs:
                if contig['id'] == contig_id:
                    length = contig['length']
            hit = {
                'contig_start': int(cols[2]),
                'contig_end': int(cols[3]),
                'plasmid_start': int(cols[4]),
                'plasmid_end': int(cols[5]),
                'plasmid': {
                    'id': cols[1],
                    'length': int(cols[6])
                },
                'coverage': float(cols[7]) / length,
                'identity': float(cols[8]) / float(cols[7])
            }
            if(hit['coverage'] >= 0.8):
                with open(log_file, "a") as fh:
                    fh.write('ref plasmids: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f' %
                        (contig_id, hit['plasmid']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity'])
                    )
                with open(args.output, "a") as fh:
                    fh.write(f"id: {contig_id} {hit}\n")
                    
                    
else:
    with open(args.output, "a") as fh:    #if not os.path.isfile(args.output):
        fh.write('None\n')

if not os.path.isfile(args.output):
    with open(args.output, "a") as fh:   
        fh.write('None\n')


        
