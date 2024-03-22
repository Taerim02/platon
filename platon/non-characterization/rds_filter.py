import functools as ft
import os
import re
import sys
import pyfastx
import argparse
import pandas as pd
import json
sys.path.append("./platon")
import platon.constants as pc

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection") 
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("tsv_file", help="TSV FASTA file")
parser.add_argument("orf_file", help="orf_file")
parser.add_argument("--mps", help="mps path")
parser.add_argument("--name", help="Output directory")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", help="Enable verbose output")

args = parser.parse_args()

def get_base_name(file_name):
    return os.path.splitext(os.path.basename(file_name))[0]
file_name = get_base_name(str(args.fasta_file))

pattern = r'_(\d+)\.fasta'
match = re.search(pattern, args.fasta_file)
match = match.group(1)

log_file = os.path.join(f'{args.name}.log')
protein_score_file = os.path.join(args.output, f'{file_name}_protein_score.txt')

contigs = {}

try:
    for record in pyfastx.Fasta(str(args.fasta_file)):
        length = len(record.seq)
        contig = {
            'id': record.name,
            'length': length,
            'sequence': str(record.seq),
            'orfs': {},
        }
        contigs[record.name] = contig
except Exception as e:
    sys.exit(f'ERROR: {str(e)}')

with open(args.orf_file, "r") as orf_file:
    orf_lines = orf_file.readlines()

for id, contig in contigs.items():
    for line in orf_lines:
        try:
            orf_dict = json.loads(line)
            if orf_dict["contig"] == id:
                orf = {
                    'start': int(orf_dict["start"]),
                    'end': int(orf_dict["end"]),
                    'strand': orf_dict["strand"],
                    'id': int(orf_dict["id"])
                }
                contig['orfs'][orf['id']] = orf
        except json.decoder.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}")
                
proteins_identified = 0
protein_id = os.path.join(args.output, 'protein_id.txt')

with open(args.tsv_file, "r") as fh:
    for line in fh:
        cols = line.split('\t')
        locus = cols[0].rpartition('_')
        contig_id = locus[0]
        orf_id = locus[2]
        if((float(cols[2]) >= pc.MIN_PROTEIN_IDENTITY) and (contig_id in contigs)):
            contig = contigs[contig_id]
            orf = contig['orfs'][int(orf_id)]
            orf['protein_id'] = cols[1]
            proteins_identified += 1
            
with open(protein_id, "a") as fh:
    fh.write(f'\tfound {proteins_identified} MPS\n')

# parse protein score file
marker_proteins = {}

with open(args.mps, "r") as fh:
    for line in fh:
        cols = line.split('\t')
        marker_proteins[cols[0]] = {
            'id': cols[0],
            'product': cols[1],
            'length': cols[2],
            'score': float(cols[3]),
        }
    
protein_score = os.path.join(args.output, f'protein_score/{file_name}_protein_score.txt')

# calculate protein score per contig
for contig in contigs.values():
    score_sum = 0.0
    for orf in contig['orfs'].values():
        if(('protein_id' in orf) and (orf['protein_id'] in marker_proteins)):
            marker_protein = marker_proteins[orf['protein_id']]
            score = marker_protein['score']
            orf['score'] = score
            orf['product'] = marker_protein['product']
            score_sum += score
    if len(contig['orfs']) > 0: 
        contig['protein_score'] = score_sum / len(contig['orfs'])
        with open(log_file, "a") as fh:
            fh.write(
                'contig RDS: contig=%s, RDS=%f, score-sum=%f, ORFs=%d\n' % 
                (contig['id'], contig['protein_score'], score_sum, len(contig['orfs']))
            )
        with open(protein_score, "a") as fh:
            fh.write(
                'contig RDS: contig=%s, RDS=%f\n' % 
                (contig['id'], contig['protein_score'])
            )
        

proteins_path = os.path.join(args.output, f'protein/{file_name}_proteins.faa')

# extract proteins from potential plasmid contigs for subsequent analyses
full_filtered_proteins_path = os.path.join(args.output, f"{args.name}_filtered.faa")

for record in pyfastx.Fasta(str(proteins_path)):
    orf_name = str(record.name).split()[0]
    contig_id = orf_name.rsplit('_', 1)[0]
    if contig_id in contigs:  # Check if contig_id is in the list of selected names
        with open(full_filtered_proteins_path, "a") as ffh:
            ffh.write(f'>{orf_name}\n')
            ffh.write(f'{record.seq}\n')


# write contig sequences to fasta files for subsequent parallel analyses
full_filtered_contig_path = os.path.join(args.output, f"{args.name}_filtered.fasta")
print(full_filtered_contig_path)
for id, contig in contigs.items():
    with open(full_filtered_contig_path, "a") as ffh:
        ffh.write(f">{contig['id']}\n")
        ffh.write(f"{contig['sequence']}\n")


    
