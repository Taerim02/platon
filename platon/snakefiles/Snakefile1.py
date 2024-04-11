import os
import glob
import re

import platon
import platon.db as db
import platon.config as cfg
import platon.functions as pf

# store variables needed in the snakefile1.py
current_path = config["current_path"]
min_protein_identity = config["min_protein_identity"]

all_file_paths = glob.glob("tmp/chunk/*.fasta")

# Get the name from the genome file
file_name = os.path.splitext(os.path.basename(cfg.genome_path))[0]

# Store the extracted file_name part
chunk_list = []
pattern = r'(.+)_(\d+)\.fasta$'

for file_path in all_file_paths:
    # Extract both the file name part and the number from the file path
    match = re.match(pattern, os.path.basename(file_path))
    if match:
        number = int(match.group(2))
        chunk_list.append(number)

rule all:
    input:
        "tmp/protein/proteins.faa", 
        "tmp/mps/protein_id.tsv",
        "tmp/orf/orf.tsv"

rule find_orfs:
    input:
        fasta_file="tmp/chunk/{file_name}_{i}.fasta", 
    
    output:
        proteins_path="tmp/protein/{file_name}_{i}_proteins.faa",
        orf_tsv = "tmp/orf/{file_name}_{i}_orf.tsv",

    params:
        verbose = cfg.verbose
    shell:
        'python {current_path}/non-characterization/orfs_find.py --verbose {params.verbose} {input.fasta_file}' 



rule count_mps:
    input:
        faa_file= "tmp/protein/{file_name}_{i}_proteins.faa"
    
    output:
        outdir = "tmp/{file_name}_{i}_diamond.tsv"
    
    params:
        threads = cfg.threads,
        tmpdir= 'tmp',
        db_path = cfg.db_path.joinpath('mps.dmnd')
    shell:
        'diamond blastp --db {params.db_path} --query {input.faa_file} --out {output.outdir} --max-target-seqs 1 --id 90 --query-cover 80 --subject-cover 80 --threads {params.threads} --tmpdir {params.tmpdir} --quiet'



rule filter_rds:
    input:
        fasta_file = "tmp/chunk/{file_name}_{i}.fasta",
        tsv_file = 'tmp/{file_name}_{i}_diamond.tsv',
        orf_tsv = "tmp/orf/{file_name}_{i}_orf.tsv",

    output:
        protein_score = "tmp/protein_score/{file_name}_{i}_protein_score.tsv",
        protein_id = 'tmp/mps/{file_name}_{i}_protein_id.tsv'
    params:
        mps = cfg.db_path.joinpath('mps.tsv'),
        outdir= 'tmp',
        verbose = cfg.verbose
    shell:         
        'python {current_path}/non-characterization/rds_filter.py --verbose {params.verbose} --mps {params.mps} --output {params.outdir} --min {min_protein_identity} {input}'
        
rule integrate_protein:
    input:
        expand('tmp/protein/{file_name}_{i}_proteins.faa', file_name=file_name, i=chunk_list) 
        
    output:
        "tmp/protein/proteins.faa"
        
    shell:
        """
        cat {input} > {output}
        """
        
rule integrate_mps:
    input:
        expand('tmp/mps/{file_name}_{i}_protein_id.tsv', file_name=file_name, i=chunk_list) 
        
    output:
        "tmp/mps/protein_id.tsv"
        
    shell:
        """
        cat {input} > {output}
        """
rule integrate_orf:
    input:
        expand("tmp/orf/{file_name}_{i}_orf.tsv", file_name=file_name, i=chunk_list)
    
    output:
        "tmp/orf/orf.tsv"
    shell:
        """
        cat {input} > {output}
        """
        