import argparse
import csv
import os 

import pyhmmer

parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("faa_file", help="Input FASTA file")
parser.add_argument("--db_path", help="db path")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()

# set a faa file and the file name  from the input file
file_name = os.path.splitext(os.path.basename(str(args.faa_file)))[0]
sequence_database_path = str(args.faa_file)  

# set the variable needed to store data processed by pyhmmer
tsv_header = ['contig', 'type', 'orf_id', 'bitscore', 'evalue']
hits_set = set()

with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
    proteins = seq_file.read_block()

# process the data with the pyhmmer
with pyhmmer.plan7.HMMFile(args.db_path) as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-100):
        for hit in hits:
            if hit.included:
                contig_id, orf_id = hit.name.decode().rsplit('_', 1)
                if f'{contig_id}_{orf_id}' not in hits_set:
                    hit_info = {
                        'type': hits.query_name.decode(),
                        'bitscore': hit.score,
                        'evalue': hit.evalue,
                    }
                    hits_set.add(f'{contig_id}_{orf_id}')
                    with open(args.output, "a", newline='') as fh:
                        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                        is_empty = fh.tell() == 0
                        if is_empty:
                            writer.writeheader()
                        tsv_row = {"contig":contig_id, "type": hit_info['type'], "bitscore":hit_info["bitscore"], "evalue":hit_info["evalue"],  "orf_id": orf_id}
                        writer.writerow(tsv_row)