import os
import pyfastx
import argparse
import csv

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection") 
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("tsv_file", help="TSV FASTA file")
parser.add_argument("orf_file", help="orf_file")
parser.add_argument("--mps", help="mps path")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help="deactivate filters; characterize all contigs")
parser.add_argument("--verbose", help="Enable verbose output")
parser.add_argument("--min", help="minimum protein identity")
args = parser.parse_args()


file_name = os.path.splitext(os.path.basename(str(args.fasta_file)))[0]

contigs = {}

for record in pyfastx.Fasta(str(args.fasta_file)):
    length = len(record.seq)
    contig = {
        'id': record.name,
        'length': length,
        'sequence': str(record.seq),
        'orfs': {}
    }
    contigs[record.name] = contig


with open(args.orf_file, "r") as orf_file:
    reader = csv.DictReader(orf_file, delimiter="\t")
    for row in reader:
        if row["contig"] in contigs:
            orf = {
                'start': int(row["start"]),
                'end': int(row["end"]),
                'strand': row["strand"],
                'id': row["id"]
            }
            contigs[row['contig']]['orfs'][orf['id']] = orf
          
proteins_identified = 0

with open(args.tsv_file, "r") as fh:
    for line in fh:
        cols = line.split('\t')
        locus = cols[0].rpartition('_')
        contig_id = locus[0]
        orf_id = locus[2]
        if((float(cols[2]) >= float(args.min)) and (contig_id in contigs)):
            contig = contigs[contig_id]
            orf = contig['orfs'][orf_id]
            orf['protein_id'] = cols[1]
            proteins_identified += 1
            
protein_id = os.path.join(args.output, f'mps/{file_name}_protein_id.tsv')
with open(protein_id, "a") as fh:
    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=['mps'])
    is_empty = fh.tell() == 0
    if is_empty:
        writer.writeheader()
    tsv_row = {'mps':proteins_identified}
    writer.writerow(tsv_row)


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

protein_score = os.path.join(args.output, f'protein_score/{file_name}_protein_score.tsv')
tsv_header = ["contig", "RDS", "score_sum", "ORFs"]


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
    contig['protein_score'] = score_sum / len(contig['orfs']) if len(contig['orfs']) > 0 else 0
    with open(protein_score, "a") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
        is_empty = fh.tell() == 0
        if is_empty:
            writer.writeheader()
        tsv_row = {"contig":contig['id'], "score_sum": score_sum, "RDS":contig['protein_score'], "ORFs":len(contig['orfs'])}
        writer.writerow(tsv_row)
