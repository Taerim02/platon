import argparse
import csv
import os 
import subprocess as sp
from pathlib import Path

import pyfastx




parser = argparse.ArgumentParser(description="Process FASTA file for ORF detection")
parser.add_argument("fasta_file", help="Input FASTA file")
parser.add_argument("--tmpdir", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument("--output", nargs='?', default=os.getcwd(), help="Output directory")
parser.add_argument('--characterize', '-c', action='store_true', help='deactivate filters; characterize all contigs')
parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
parser.add_argument("--min", help="minimum circular basepair overlap")

args = parser.parse_args()

file_name = os.path.splitext(os.path.basename(str(args.fasta_file)))[0]

contigs = {}

for record in pyfastx.Fasta(str(args.fasta_file)):
    length = len(record.seq)
    contig = {
        'id': record.name,
        'length': length,
        'sequence': str(record.seq),
    }
    contigs[record.name] = contig

tmpdir_path = Path(args.tmpdir)

def write_sequence_to_file(contig, contig_split_position, header):
    file_path = tmpdir_path.joinpath(f"{contig['id']}-{header}.fasta")
    if header == 'a':
        sequence = contig['sequence'][:contig_split_position]
    elif header == 'b':
        sequence = contig['sequence'][contig_split_position:]
    with file_path.open(mode='w') as fh:
        fh.write('>{header}\n')
        fh.write(sequence +'\n ')
    return file_path, sequence

tsv_header = ['contig', 'length', 'mismatches', 'prime5End', 'prime3Start', 'prime3End', 'seq-a-len', 'seq-b-len']

for id, contig in contigs.items():
    contig_split_position = int(contig['length'] / 2)
    contig_fragment_a_path, contig_fragment_a_seq = write_sequence_to_file(contig, contig_split_position, 'a')
    contig_fragment_b_path, contig_fragment_b_seq = write_sequence_to_file(contig, contig_split_position, 'b')
    cmd = [
        'nucmer',
        '-f',  # only forward strand
        '-l', '40',  # increase min match length to 40 bp
        '--threads=1',
        '-p', contig['id'],
        str(contig_fragment_b_path),
        str(contig_fragment_a_path)
    ]
    proc = sp.run(
        cmd,
        cwd=str(tmpdir_path),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )

    if(proc.returncode != 0):
        print('failed')
        
    has_match = False
    with tmpdir_path.joinpath(f"{contig['id']}.delta").open() as fh:
        for line in fh:
            line = line.rstrip()
            if(line[0] == '>'):
                has_match = True
            elif(has_match):
                cols = line.split(' ')
                if(len(cols) == 7):
                    start_b = int(cols[0])
                    end_b = int(cols[1])
                    start_a = int(cols[2])
                    end_a = int(cols[3])
                    mismatches = int(cols[4])
                    alignment_a = end_a - start_a + 1
                    alignment_b = end_b - start_b + 1
                    if(alignment_a == alignment_b
                        and alignment_a > int(args.min)
                        and (mismatches / alignment_a) < 0.05
                        and end_b == len(contig_fragment_b_seq)
                        and start_a == 1):
                        link = {
                            'length': alignment_a,
                            'mismatches': mismatches,
                            'prime5End': alignment_a,
                            'prime3Start': contig['length'] - alignment_b + 1,
                            'prime3End': contig['length'],
                        }
                        with open(args.output, "a", newline='') as fh:
                            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=tsv_header)
                            is_empty = fh.tell() == 0
                            if is_empty:
                                writer.writeheader()
                            tsv_row = {'contig':contig['id'], 'length':link['length'], 'mismatches':link['mismatches'],
                                'prime5End':link['prime5End'], 'prime3Start':link['prime3Start'], 'prime3End':link['prime3End'], 
                                'seq-a-len':len(contig_fragment_a_seq), 'seq-b-len':len(contig_fragment_b_seq)}
                            writer.writerow(tsv_row)
                        break

