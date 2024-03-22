import os
import pyfastx
import argparse
from pathlib import Path
import subprocess as sp
import sys

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("contig_file", help="Input FASTA file")
parser.add_argument("rrnas_output", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()

log_file = os.path.join(f'{args.name}.log')

contig_plasmid_hits = []


if not os.path.exists('tmp/rrnas'):
    os.mkdir('tmp/rrnas')
    
if os.path.getsize(args.rrnas_output) > 0:
    with open(args.rrnas_output, "r") as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                print(cols)
                hit = {
                    'type': cols[0],
                    'start': int(cols[7]),
                    'end': int(cols[8]),
                    'strand': cols[9],
                    'bitscore': float(cols[14]),
                    'evalue': float(cols[15])
                }
                contig_plasmid_hits.append(hit)
                with open(log_file, "a") as fh:
                    fh.write(
                    'rRNAs: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s' %
                    (cols[2], hit['type'], hit['start'], hit['end'], hit['strand'])
                )
                with open(args.output, "a") as fh:  ## how to address the contig_id?
                    fh.write(f"id: {cols[2]} {hit}\n")                          
            if len(contig_plasmid_hits) > 0:
                with open(log_file, "a") as fh:
                    fh.write('rRNAs: contig=%s, rRNAs=%s' % ('id', len(contig_plasmid_hits)))

else:
    with open(args.output, "a") as fh:
          fh.write('None\n')


if not os.path.isfile(args.output):
    with open(args.output, "a") as fh:   
        fh.write('None\n')





