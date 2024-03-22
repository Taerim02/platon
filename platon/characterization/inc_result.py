import os, sys
import pyfastx
from pathlib import Path
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(description="Process FASTA file for Incompatibility groups detection")
parser.add_argument("contig_file", help="Input FASTA file")
parser.add_argument("inc_output", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()


log_file = os.path.join(f'{args.name}.log')

hits_per_pos = {}
contig_id = ''

if os.path.getsize(args.inc_output) > 0:
    with open(args.inc_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            new_contig_id = cols[1]
            if contig_id != new_contig_id:
                if len(list(hits_per_pos.values())) > 0:
                    with open(log_file, "a") as fh:
                        fh.write('inc-type: contig=%s, inc-types=%s\n'  % (contig_id, len(list(hits_per_pos.values()))))
                hits_per_pos = {}
                contig_id = new_contig_id
            hit = {
                'type': cols[0],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'strand': '+' if cols[4] == 'plus' else '-',
                'identity': float(cols[5]) / 100,
                'coverage': float(cols[6]) / 100,
                'bitscore': int(cols[7])
            }
            if(hit['coverage'] >= 0.6):
                hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                if(hit_pos in hits_per_pos):
                    former_hit = hits_per_pos[hit_pos]
                    if(hit['bitscore'] > former_hit['bitscore']):
                        hits_per_pos[hit_pos] = hit
                        with open(log_file, "a") as fh:
                            fh.write(
                                'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s\n' %
                                (contig_id, hit['type'], hit['start'], hit['end'], hit['strand'])
                            )
                        with open(args.output, "a") as fh:
                            fh.write(f"id: {contig_id} {hit}\n")
                else:
                    hits_per_pos[hit_pos] = hit
                    with open(log_file, "a") as fh:
                        fh.write(
                            'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s\n' %
                            (contig_id, hit['type'], hit['start'], hit['end'], hit['strand'])
                        )
                    with open(args.output, "a") as fh:
                        fh.write(f"id: {contig_id} {hit}\n")        

else: 
    with open(args.output, "a") as fh:
        fh.write('None\n')
        
if not os.path.isfile(args.output):
    with open(args.output, "a") as fh:   
        fh.write('None\n')
