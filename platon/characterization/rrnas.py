import argparse
import csv
import os

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("rrnas_output", help="Input FASTA file")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()

# set the variable needed to store data
tsv_header = ['contig', 'type', 'start', 'end', 'strand', 'bitscore', 'evalue']

if os.path.getsize(args.rrnas_output) > 0:
    with open(args.rrnas_output, "r") as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                hit_info = {
                    'type': cols[0],
                    'start': int(cols[7]),
                    'end': int(cols[8]),
                    'strand': cols[9],
                    'bitscore': float(cols[14]),
                    'evalue': float(cols[15])
                }
                with open(args.output, "a", newline='') as fh:
                    writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                    is_empty = fh.tell() == 0
                    if is_empty:
                        writer.writeheader()
                    tsv_row = {"contig":cols[2], "type": hit_info['type'], "start":hit_info["start"], "end":hit_info["end"], 
                        "strand":hit_info["strand"], "bitscore": hit_info["bitscore"], "evalue":hit_info["evalue"]}
                    writer.writerow(tsv_row)
