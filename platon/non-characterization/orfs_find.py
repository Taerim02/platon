import functools as ft
import os
import re
import sys
import pyfastx
import pyrodigal
import argparse
import pandas as pd
import json


parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument("--verbose", help="Enable verbose output")
args = parser.parse_args()

file_name = os.path.splitext(os.path.basename(str(args.fasta_file)))[0]

pattern = r'_(\d+)\.fasta'
match = re.search(pattern, args.fasta_file)
match = match.group(1)

contigs = {}
raw_contigs = []
try:
    for record in pyfastx.Fasta(str(args.fasta_file)):
        length = len(record.seq)
        contig = {
            'id': record.name,
            'length': length,
            'sequence': str(record.seq),
            'orfs': {},
            'is_circular': False,
            'inc_types': [],
            'amr_hits': [],
            'mobilization_hits': [],
            'orit_hits': [],
            'replication_hits': [],
            'conjugation_hits': [],
            'rrnas': [],
            'plasmid_hits': [],
            'type': []
        }
        raw_contigs.append(contig)
        contigs[record.name] = contig
        
except Exception as e:
    sys.exit(f'ERROR: {str(e)}')
    

orf_txt =  f'tmp/orf/{file_name}_orf.txt'
proteins_path = f'tmp/protein/{file_name}_proteins.faa'
log_file = os.path.join(f'{args.name}.log')



for record in pyfastx.Fasta(str(args.fasta_file)):
    orf_finder = pyrodigal.GeneFinder(meta=True, closed=True, mask=True)
    genes = orf_finder.find_genes(str(record.seq))
    with open(proteins_path, "a") as dst:
        genes.write_translations(dst, sequence_id=record.name)
    for i, gene in enumerate(genes):
        orf_id = i + 1
        strand_match = re.search(r"strand=([-+])", str(gene))
        if strand_match: 
            orf = {
                'start': gene.begin,
                'end': gene.end,
                'strand': strand_match.group(1),
                'id': orf_id
            }
            contig = contigs.get(record.name, None)
            if contig is not None:
                contig['orfs'][orf_id] = orf
                with open(log_file, "a") as fh:  # Changed to append mode ("a")
                    fh.write(
                        f'ORFs: found! contig={contig["id"]}, start={orf["start"]}, end={orf["end"]}, strand={orf["strand"]}\n'
                    )
                with open(orf_txt, "a") as fh:
                    fh.write(f'{{"contig":"{contig["id"]}", "start":{orf["start"]}, "end":{orf["end"]}, "strand":"{orf["strand"]}", "id":{orf["id"]}}}\n')

if not os.path.exists(proteins_path):  # Check if proteins_path file was created
    sys.exit('Error: ORF prediction failed!')

no_orfs = sum(len(contigs[k]['orfs']) for k in contigs)  # Calculate total number of ORFs

if(args.verbose):
    print(f'\tfound {no_orfs} ORFs')

with open(log_file, "a") as fh:  # Changed to append mode ("a")
    fh.write(f'ORF detection: # ORFs={no_orfs}\n')

