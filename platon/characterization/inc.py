import argparse
import csv
import os

"""Search for incompatibility motifs in the metagenomic module of platon."""

parser = argparse.ArgumentParser(description="Process FASTA Incompatibility groups detection")
parser.add_argument("inc_output", help="Input FASTA file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")

args = parser.parse_args()


# set the variable needed to store data
tsv_header = ['contig', 'type', 'start', 'end', 'strand', 'identity', 'coverage', 'bitscore', 'hit_pos']
hits_per_pos = {}

# extract information and store into tsv file
if os.path.getsize(args.inc_output) > 0:
    with open(args.inc_output,"r") as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            contig_id = cols[1]
            hit = {
                'type': cols[0],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'strand': '+' if cols[4] == 'plus' else '-',
                'identity': float(cols[5]) / 100,
                'coverage': float(cols[6]) / 100,
                'bitscore': float(cols[7])
            }
            if(hit['coverage'] >= 0.6):
                hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                hit['hit_pos'] = hit_pos
                with open(args.output, "a", newline='') as fh:
                    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                    is_empty = fh.tell() == 0
                    if is_empty:
                        writer.writeheader()
                    if(hit_pos in hits_per_pos):
                        former_hit = hits_per_pos[hit_pos]
                        if(hit['bitscore'] > former_hit['bitscore']):
                            tsv_row = {"contig":contig_id, "type": hit['type'], "start":hit["start"], "end": hit["end"], 'strand': hit['strand'],
                                'identity': hit['identity'], 'coverage':hit['coverage'], "bitscore":hit["bitscore"], "hit_pos":hit["hit_pos"]}
                            writer.writerow(tsv_row)
                    else:
                        tsv_row = {"contig":contig_id, "type": hit['type'], "start":hit["start"], "end": hit["end"], 'strand': hit['strand'],
                                'identity': hit['identity'], 'coverage':hit['coverage'], "bitscore":hit["bitscore"], "hit_pos":hit["hit_pos"]}
                        writer.writerow(tsv_row)
            
                        


        



