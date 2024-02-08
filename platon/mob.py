import pyhmmer
import os 
import argparse


parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("faa_file", help="Input FASTA file")
parser.add_argument("--name", help="original fasta file name for generating a log file")
parser.add_argument("--db_path", help="db path")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

args = parser.parse_args()


log_file = os.path.join(f'{args.name}.log')
file_name = os.path.splitext(os.path.basename(str(args.faa_file)))[0]
sequence_database_path = str(args.faa_file)  

contig_dict = {}
hits_set = set()

with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
    proteins = seq_file.read_block()

with pyhmmer.plan7.HMMFile(args.db_path) as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-10):
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
                    if contig_id in contig_dict:
                        contig_dict[contig_id].append(hit_info)
                    else:
                        contig_dict[contig_id] = [hit_info]
                    with open(log_file, "a") as fh:
                        fh.write('mob genes: hit! contig=%s, type=%s' %
                            (contig_id, hit_info['type'])
                        )
                    with open(args.output, "a") as fh:
                            fh.write(f"id: {contig_id} orf_id: {orf_id} {hit_info}\n")


with open(log_file, "a") as fh:
    for contig_id in contig_dict.keys():
        fh.write(f'mob genes: contig={contig_id}, mobs={len(contig_dict[contig_id])}\n')
    

if not os.path.isfile(args.output):
    with open(args.output, "a") as fh:
        fh.write(f'None\n') 
