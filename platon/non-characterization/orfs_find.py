import os, re, sys
import pyfastx
import pyrodigal
import argparse
import csv

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument("--verbose", help="Enable verbose output")


args = parser.parse_args()

file_name = os.path.splitext(os.path.basename(str(args.fasta_file)))[0]

contigs = {}
try:
    for record in pyfastx.Fasta(str(args.fasta_file)):
        length = len(record.seq)
        contig = {
            'id': record.name,
            'length': length,
            'sequence': str(record.seq),
            'orfs': {}
        }
        contigs[record.name] = contig
        
except Exception as e:
    sys.exit(f'ERROR: {str(e)}')
    
proteins_path = f'tmp/protein/{file_name}_proteins.faa'
orf_tsv =  f'tmp/orf/{file_name}_orf.tsv'
tsv_header = ["contig", "start", "end", "strand", "id"]

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
                with open(orf_tsv, "a", newline='') as fh:
                    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                    is_empty = fh.tell() == 0
                    if is_empty:
                        writer.writeheader()
                    tsv_row = {"contig":contig["id"], "start":orf["start"], "end":orf["end"], "strand":orf["strand"], "id":orf["id"]}
                    writer.writerow(tsv_row)
                    

# Check if proteins_path file was created
if not os.path.exists(proteins_path):  
    sys.exit('Error: ORF prediction failed!')



