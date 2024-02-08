import functools as ft
import json
import logging
import os, sys, re
from pathlib import Path
import pyfastx
import shutil
from snakemake import snakemake
import pandas as pd
import ast


import platon
import platon.__init__ as init
import platon.db as db
import platon.config as cfg
import platon.constants as pc
import platon.functions as pf
import platon.utils as pu


def main(raw_contigs, contigs, filtered_contigs, args, log, output_path):
        
    if not os.path.exists(output_path.joinpath('tmp')):    
            os.mkdir(output_path.joinpath('tmp'))
            
    for folder in [ "rrnas", "orit", "inc", "ref", 'orf', 'function', 'chunk', 'protein']:
        if not os.path.exists(output_path.joinpath(f'tmp/{folder}')):
            os.mkdir(output_path.joinpath(f'tmp/{folder}'))
            
    
    
    for folder in ["amr", "rep", "mob", "conj", "cir", "rrnas", "orit", "inc", "ref"]:
        if not os.path.exists(output_path.joinpath(f'tmp/function/{folder}')):
            os.mkdir(output_path.joinpath(f'tmp/function/{folder}'))
    
    pf.contigs_into_chunks(contigs, pc.DEFAULT_CONTIG_SIZE, output_path) 
    
    result1 = snakemake(
        snakefile="platon/Snakefile1",
        workdir=output_path,
        config={"current_path": Path.cwd()},
        cores=cfg.threads, 
    )
    
    
    name = pf.get_base_name(cfg.genome_path)
    
    full_filtered_contig_path = os.path.join(output_path.joinpath('tmp'), f"{name}_filtered.fasta")
    pf.fasta_into_chunk(full_filtered_contig_path, pc.DEFAULT_CONTIG_SIZE, output_path)
    
    
    full_filtered_proteins_path = os.path.join(output_path.joinpath('tmp'), f"{name}_filtered.faa")
    pf.faa_into_chunk(full_filtered_proteins_path, pc.DEFAULT_CONTIG_SIZE, output_path)
    
    
    
    result2 = snakemake(
        snakefile="platon/Snakefile2",
        workdir=output_path,
        config={"current_path": Path.cwd()},
        cores=cfg.threads, scheduler="greedy",
    )
    
    
    orf_files = os.listdir(output_path.joinpath('tmp/orf'))
    tmp_contigs = {}
    for file in orf_files:
        with open(os.path.join(output_path.joinpath('tmp/orf'), file),"r") as fh:
            for line in fh:
                line = json.loads(line) 
                contig_id = line['contig']
                orf_id = line['id']
                contigs[contig_id]['orfs'][orf_id]= line
                tmp_contigs[contig_id] = contigs[contig_id]
    contigs = tmp_contigs
    
    protein_score_pattern = r'contig RDS: contig=([a-zA-Z0-9_.]+), RDS=(-?[0-9.]+)'
    protein_score_files = os.listdir(output_path.joinpath('tmp/protein_score'))
    for f in protein_score_files:
        with open(os.path.join(output_path.joinpath('tmp/protein_score'), f),"r") as fh:
            for line in fh:
                match = re.match(protein_score_pattern, line)
                if match:
                    contig_id = match.group(1)
                    if contig_id in contigs:
                        contigs[contig_id]['protein_score']= float(match.group(2))
    
                    
    hmm = ['amr_hits', 'conjugation_hits', 'mobilization_hits', 'replication_hits']                  

    for names in [ ['amr.txt','amr_hits'], ['rrnas.txt','rrnas'], ['orit.txt','orit_hits'], ['cir.txt','is_circular'], ['inc.txt','inc_types'], 
            ['ref.txt','plasmid_hits'], ['mob.txt','mobilization_hits'], ['rep.txt','replication_hits'], ['conj.txt','conjugation_hits']]:
            
        if names[1] in hmm: 
            pf.extract_function_info_hmm(contigs, names[0], names[1], output_path)
        else:
            pf.extract_function_info(contigs, names[0], names[1], output_path)
         
    if(args.characterize):  # skip protein score based filtering
        filtered_contigs = contigs
    elif(args.mode == 'sensitivity'):  # skip protein score based filtering but apply rRNA filter
        filtered_contigs = {k: v for (k, v) in contigs.items() if v['rrnas'] == 0}
    elif(args.mode == 'specificity'):
        filtered_contigs = {k: v for (k, v) in contigs.items() if v['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD}
    else:
        filtered_contigs = {k: v for (k, v) in contigs.items() if pf.filter_contig(v)}
        
        
    # lookup AMR genes
    amr_genes = {}
    with cfg.db_path.joinpath('ncbifam-amr.tsv').open() as fh:
        for line in fh:
            cols = line.split('\t')
            amr_genes[cols[0]] = {
                'gene': cols[4],
                'product': cols[8]
            }
    
    for id, contig in contigs.items():
        for hit in contig['amr_hits']:
            amr_gene = amr_genes[hit['hmm-id']]
            hit['gene'] = amr_gene['gene']
            hit['product'] = amr_gene['product']
    
    
    # print results to tsv file and STDOUT
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.tsv')
    log.debug('output: tsv=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        fh.write(pc.HEADER + '\n')
        if(len(filtered_contigs) > 0):
            print(pc.HEADER)
            for id in sorted(filtered_contigs, key=lambda k: -filtered_contigs[k]['length']):
                c = filtered_contigs[id]
                cov = 'NA' if c['coverage'] == 0 else f"{c['coverage']:4.1f}"
                circ = 'yes' if c['is_circular'] else 'no'
                line = f"{c['id']}\t{c['length']}\t{cov}\t{len(c['orfs'])}\t{c['protein_score']:3.1f}\t{circ}\t{len(c['inc_types'])}\t{len(c['replication_hits'])}\t{len(c['mobilization_hits'])}\t{len(c['orit_hits'])}\t{len(c['conjugation_hits'])}\t{len(c['amr_hits'])}\t{len(c['rrnas'])}\t{len(c['plasmid_hits'])}"
                print(line)
                log.info(
                    'plasmid: id=%s, len=%d, cov=%s, ORFs=%d, RDS=%f, circ=%s, incs=%s, rep=%d, mob=%d, oriT=%d, con=%d, AMRs=%d, rRNAs=%d, refs=%d' %
                    (c['id'], c['length'], cov, len(c['orfs']), c['protein_score'], c['is_circular'], len(c['inc_types']), len(c['replication_hits']), len(c['mobilization_hits']), len(c['orit_hits']), len(c['conjugation_hits']), len(c['amr_hits']), len(c['rrnas']), len(c['plasmid_hits']))
                )
                fh.write(f'{line}\n')
        else:
            print('No potential plasmid contigs found!')
            print(pc.HEADER)
                  
    shutil.rmtree(output_path.joinpath('tmp'))
    return raw_contigs, contigs, filtered_contigs
