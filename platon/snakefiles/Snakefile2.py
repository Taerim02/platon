import os, re 
import glob
import pyfastx

import platon
import platon.db as db
import platon.config as cfg
import platon.functions as pf

file_name = pf.get_base_name(cfg.genome_path)
current_path = config["current_path"]
min_circ_basepair_overlap = config["min_circ_basepair_overlap"]

# Define the regular expression pattern
pattern_faa = r'(.+)_(\d+)_filtered.faa$'
pattern_fasta = r'(.+)_(\d+)_filtered.fasta$'
all_fa_file_paths = glob.glob('tmp/*_filtered.faa')
all_fasta_file_paths = glob.glob('tmp/*_filtered.fasta')


# Store the extracted file_name part
faa_chunk_list = []
faa_files = []
for file_path in all_fa_file_paths:
    # Extract both the file name part and the number from the file path
    match = re.match(pattern_faa, os.path.basename(file_path))   
    if match:
        faa_files.append(file_path)
        number = int(match.group(2))
        faa_chunk_list.append(number)

fasta_chunk_list = []
fasta_files = []

for file_path in all_fasta_file_paths:
    # Extract both the file name part and the number from the file path
    match = re.match(pattern_fasta, os.path.basename(file_path))   
    if match:
        fasta_files.append(file_path)
        number = int(match.group(2))
        fasta_chunk_list.append(number)        


def get_blast_word_size(contig_name):
    config_file = f"tmp/contigs/{contig_name}.fasta"
    fasta = pyfastx.Fasta(config_file)
    record = fasta[0]
    length = len(record.seq)
    blast_word_size = int(length / 1000)
    if(blast_word_size < 11):
        blast_word_size = 11
    return blast_word_size


current_directory = os.getcwd()
all_contigs = glob.glob(os.path.join(current_directory, 'tmp/contigs/*.fasta'))
contig_names = []

for contig_file in all_contigs:
    # Extract both the file name part and the number from the file path
    match = re.match(r'(.+)\.fasta$', os.path.basename(contig_file))   
    if match:
        contig_names.append(match.group(1))

faa = ["amr", "rep", "mob", "conj"]  # List of possible values for the tsv wildcard
fasta = ["cir", "rrnas", "orit", "inc", "ref"]


rule all:
    input:
        'tmp/function/amr.tsv', 'tmp/function/rep.tsv' , 'tmp/function/mob.tsv' ,'tmp/function/conj.tsv' ,'tmp/function/ref.tsv' ,
        'tmp/function/cir.tsv' ,'tmp/function/inc.tsv' ,'tmp/function/orit.tsv' , 'tmp/function/rrnas.tsv'


rule search_amr_genes:
    input:
        faa_file = 'tmp/{file_name}_{i}_filtered.faa',
    output:
        amr = touch('tmp/function/amr/{file_name}_{i}_amr.tsv')   
    params:
        db_path = cfg.db_path.joinpath('ncbifam-amr'),
        name = file_name
    shell:
        'python {current_path}/characterization/amr.py --db_path {params.db_path} --name {params.name} --output {output.amr} {input.faa_file}'


rule search_replication_genes: 
    input:
        faa_file = 'tmp/{file_name}_{i}_filtered.faa',
    output:
        rep = touch('tmp/function/rep/{file_name}_{i}_rep.tsv')  
    params:
        db_path = cfg.db_path.joinpath('replication')
    shell:
        'python {current_path}/characterization/rep.py --db_path {params.db_path} --output {output.rep} {input.faa_file}'


rule search_mobilization_genes:
    input:
        faa_file = 'tmp/{file_name}_{i}_filtered.faa',
    output:
        mob = touch('tmp/function/mob/{file_name}_{i}_mob.tsv')  
    params:
        db_path = cfg.db_path.joinpath('mobilization')
    shell:
        'python {current_path}/characterization/mob.py --db_path {params.db_path} --output {output.mob} {input.faa_file}'

rule search_conjugation_genes:
    input:
        faa_file = 'tmp/{file_name}_{i}_filtered.faa',
    output:
        conj = touch('tmp/function/conj/{file_name}_{i}_conj.tsv')   
    params:
        db_path = cfg.db_path.joinpath('conjugation'),
    shell:
        'python {current_path}/characterization/conj.py --db_path {params.db_path} --output {output.conj} {input.faa_file}'


rule test_circularity:
    input:
        fasta_file = "tmp/{file_name}_{i}_filtered.fasta"
    output:
        cir = touch('tmp/function/cir/{file_name}_{i}_cir.tsv') 
    params:
        tmpdir = str(cfg.tmp_path),
    shell: 
        'python {current_path}/characterization/circularity.py --tmpdir {params.tmpdir} --output {output.cir} --min {min_circ_basepair_overlap} {input.fasta_file}'



rule search_rrnas:
    input:
        fasta_file = "tmp/{file_name}_{i}_filtered.fasta"
    output:
        rrnas = touch('tmp/function/rrnas/{file_name}_{i}_rrnas.tsv')
    params:
        db_path = str(cfg.db_path.joinpath('rRNA')),
        rrnas = "tmp/rrnas/{file_name}_{i}.rrnas.blast.out",
    shell: 
        r"""
        cmscan \
        --noali \
        --cut_tc \
        --cpu 1 \
        --tblout {params.rrnas} \
        {params.db_path} \
        {input.fasta_file} > tmp/rrnas/output.txt

        python {current_path}/characterization/rrnas_result.py --output {output.rrnas} {params.rrnas} 
        """
        
rule search_inc_type:
    input:
        contig_file = "tmp/{file_name}_{i}_filtered.fasta"
    output:
        inc = touch('tmp/function/inc/{file_name}_{i}_inc.tsv')
    params:
        db_path = str(cfg.db_path.joinpath('inc-types.fasta')),
        inc = "tmp/inc/{file_name}_{i}.inc.blast.out",
        name = file_name
    shell:
        r"""
        blastn \
        -query {params.db_path} \
        -subject {input.contig_file} \
        -num_threads 1 \
        -culling_limit 1 \
        -perc_identity 90 \
        -outfmt '6 qseqid sseqid sstart send sstrand pident qcovs bitscore' \
        -out {params.inc}

        python {current_path}/characterization/inc_result.py --name {params.name} --output {output.inc} {params.inc}
        """

rule search_orit_sequences:
    input:
        contig_file = "tmp/{file_name}_{i}_filtered.fasta"
    output:
        orit = touch('tmp/function/orit/{file_name}_{i}_orit.tsv')
    params:
        db_path = str(cfg.db_path.joinpath('orit')),
        orit = "tmp/orit/{file_name}_{i}.orit.blast.out",
    shell:
        r"""
        blastn \
        -query {input.contig_file} \
        -db {params.db_path} \
        -num_threads 1 \
        -culling_limit 1 \
        -perc_identity 90 \
        -evalue 1E-5 \
        -outfmt '6 qseqid sseqid qstart qend sstart send slen length nident' \
        -out {params.orit}

        python {current_path}/characterization/orit_result.py --output {output.orit} {params.orit}
        """

rule search_reference_plasmids:  
    input:
        contig_file = "tmp/{file_name}_{i}_filtered.fasta"
    output:
        ref = touch('tmp/function/ref/{file_name}_{i}_ref.tsv')
    params:
        db_path = str(cfg.db_path.joinpath('refseq-plasmids')),
        ref = "tmp/ref/{file_name}_{i}.refplas.blast.out",  #{params.word_size}
    shell:
        r"""
        blastn \
        -query {input.contig_file} \
        -db {params.db_path} \
        -num_threads 1 \
        -culling_limit 1 \
        -perc_identity 80 \
        -outfmt '6 qseqid sseqid qstart qend sstart send slen length nident' \
        -out {params.ref}

        python {current_path}/characterization/ref_result.py --output {output.ref} {input.contig_file} {params.ref}
        """

rule integrate_amr_files:
    input:
       expand("tmp/function/amr/{file_name}_{i}_amr.tsv", file_name=file_name, i=faa_chunk_list)  
    output:
        'tmp/function/amr.tsv'   
    shell:
        """
        cat {input} > {output}
        """

rule integrate_mob_files:
    input:
       expand("tmp/function/mob/{file_name}_{i}_mob.tsv", file_name=file_name, i=faa_chunk_list)
    output:
        'tmp/function/mob.tsv'   
    shell:
        """
        cat {input} > {output}
        """

rule integrate_conj_files:
    input:
       expand("tmp/function/conj/{file_name}_{i}_conj.tsv", file_name=file_name, i=faa_chunk_list) 
    output:
        'tmp/function/conj.tsv'  
    shell:
        """
        cat {input} > {output}
        """

rule integrate_cir_files:
    input:
       expand("tmp/function/cir/{file_name}_{i}_cir.tsv", file_name=file_name, i=fasta_chunk_list)  
    output:
        'tmp/function/cir.tsv'   
    shell:
        """
        cat {input} > {output}
        """

rule integrate_rrnas_files:
    input:
       expand("tmp/function/rrnas/{file_name}_{i}_rrnas.tsv", file_name=file_name, i=fasta_chunk_list)  
    output:
        touch('tmp/function/rrnas.tsv')   
    shell:
        """
        cat {input} > {output}
        """        

rule integrate_orit_files:
    input:
       expand("tmp/function/orit/{file_name}_{i}_orit.tsv", file_name=file_name, i=fasta_chunk_list)   
    output:
        'tmp/function/orit.tsv'  
    shell:
        """
        cat {input} > {output}
        """        
    
rule integrate_ref_file:
    input:
       expand("tmp/function/ref/{file_name}_{i}_ref.tsv", file_name=file_name, i=fasta_chunk_list)  
    output:
        'tmp/function/ref.tsv'   
    shell:
        """
        cat {input} > {output}
        """       

rule integrate_rep_files:
    input:
       expand("tmp/function/rep/{file_name}_{i}_rep.tsv", file_name=file_name, i=faa_chunk_list)  
    output:
        'tmp/function/rep.tsv'   
    shell:
        """
        cat {input} > {output}
        """        

rule integrate_inc_files:
    input:
       expand("tmp/function/inc/{file_name}_{i}_inc.tsv", file_name=file_name, i=fasta_chunk_list)  
    output:
        'tmp/function/inc.tsv'   
    shell:
        """
        cat {input} > {output}
        """            

