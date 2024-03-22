import os
import pyfastx
import argparse
from pathlib import Path
import subprocess as sp
import sys

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("contig_file", help="Input FASTA file")
parser.add_argument("orit_output", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()


log_file = os.path.join(f'{args.name}.log')


contig_id = ''
contig_plasmid_hits = []

if os.path.getsize(args.orit_output) > 0:
    with open(args.orit_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            new_contig_id = cols[0]
            if contig_id != new_contig_id:
                if len(contig_plasmid_hits) > 0:
                    with open(log_file, "a") as fh:
                        fh.write('oriT: contig=%s, oriT=%s' % (contig_id, len(contig_plasmid_hits)))   
                contig_plasmid_hits = []
                contig_id = new_contig_id
            hit = {
                'contig_start': int(cols[2]),
                'contig_end': int(cols[3]),
                'orit_start': int(cols[4]),
                'orit_end': int(cols[5]),
                'orit': {
                    'id': cols[1],
                    'length': int(cols[6])
                },
                'coverage': float(cols[7]) / int(cols[6]),
                'identity': float(cols[8]) / float(cols[7])
            }
            if(hit['coverage'] >= 0.9 and hit['identity'] >= 0.9):
                contig_plasmid_hits.append(hit)
                with open(log_file, "a") as fh:
                    fh.write(
                        'oriT: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f' %
                        (contig_id, hit['orit']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity'])
                    )
                with open(args.output, "a") as fh:
                    fh.write(f"id: {contig_id} {hit}\n")  
                    
else:
    with open(args.output, "a") as fh:
        fh.write('None\n')
        
if not os.path.isfile(args.output):
    with open(args.output, "a") as fh:   
        fh.write('None\n')

