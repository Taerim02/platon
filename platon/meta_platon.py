import os
from pathlib import Path
import shutil
from snakemake import snakemake
import platon
import platon.config as cfg
import platon.constants as pc
import platon.functions as pf
import csv

def main(raw_contigs, contigs, args, log, output_path):
    '''
    This is the metagenomic module of platon to scale better for large-size datasets.
    Snakemake is used to help parallel computing in order to improve running time.
    '''
        
    # prepare folders to gather results for metagenomic module.
    if not os.path.exists(output_path.joinpath('tmp')):    
        os.mkdir(output_path.joinpath('tmp'))

        
    for folder in [ "rrnas", "orit", "inc", "ref", 'orf', 'function', 'chunk', 'protein', 'protein_score', 'mps', 'log']:
        if not os.path.exists(output_path.joinpath(f'tmp/{folder}')):
            os.mkdir(output_path.joinpath(f'tmp/{folder}'))

            
    for folder in ["amr", "rep", "mob", "conj", "cir", "rrnas", "orit", "inc", "ref"]:
        if not os.path.exists(output_path.joinpath(f'tmp/function/{folder}')):
            os.mkdir(output_path.joinpath(f'tmp/function/{folder}'))


    # divide a fasta file into chunks according to the default contig size.
    pf.contigs_into_chunks(contigs, pc.DEFAULT_CONTIG_SIZE, output_path)


    # first Snakemake session starts for non-characsteristic part. 
    result1 = snakemake.snakemake(
        snakefile="snakefiles/Snakefile1",
        workdir= output_path,
        config={"current_path": Path.cwd(), "min_protein_identity": pc.MIN_PROTEIN_IDENTITY},
        cores=cfg.threads, quiet=True
    )

    if result1:
        if(args.verbose):
                print("Workflow for non-characterization completed successfully!")
        log.info("Workflow for non-characterization completed.")
    else:
        if(args.verbose):
            print("Workflow for non-characterization failed.")
        log.debug("Workflow for non-characterization failed.")

    # predict open reading frames with Pyrodigal.
    if(args.verbose):
      print('predict ORFs...')

    # parse the result from pyrodigal: ORF.
    with open(output_path.joinpath('tmp/orf/orf.tsv'), "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["contig"] in contigs:
                orf = {
                    'start': int(row["start"]),
                    'end': int(row["end"]),
                    'strand': row["strand"],
                    'id': row["id"]
                }
                contigs[row['contig']]['orfs'][int(row['id'])] = orf
                
    # Calculate total number of ORFs
    no_orfs = sum(len(contigs[k]['orfs']) for k in contigs)  
    
    if(args.verbose):
        print(f'\tfound {no_orfs} ORFs')
    
    log.info('ORF detection: # ORFs=%d', no_orfs)
    
    # retain / exclude contigs without ORFs
    if(args.characterize or args.mode == 'sensitivity' or args.mode == 'accuracy'):
        log.info('ORF contig filter disabled! # passed contigs=%s', len(contigs))
    else:  # exclude contigs without ORFs in specificity mode
        tmp_contigs = {}
        for contig in filter(lambda k: len(k['orfs']) > 0, contigs.values()):
            tmp_contigs[contig['id']] = contig
        no_removed_contigs = len(contigs) - len(tmp_contigs)
        if(args.verbose):
            print(f'\tdiscarded {no_removed_contigs} contig(s) without ORFs')
        contigs = tmp_contigs
        log.info('ORF contig filter: # discarded=%s, # remaining=%s', no_removed_contigs, len(contigs))
        
    
    if(args.verbose):
      print('search marker protein sequences (MPS)...')
      
    # Calculate total number of MPS
    no_mps = 0
    with open(output_path.joinpath("tmp/mps/protein_id.tsv"),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            no_mps += int(row["mps"])
    
    if(args.verbose):
        print(f'\tfound {no_mps} MPS')
    
    log.info('MPS detection: # MPS=%d', no_mps)
    
    if(args.verbose):
        print('compute replicon distribution scores (RDS)...')
                   
    protein_score_files = os.listdir(output_path.joinpath('tmp/protein_score'))
    for fh in protein_score_files:
        with open(os.path.join(output_path.joinpath('tmp/protein_score'), fh),"r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig_id = row["contig"]
                contigs[contig_id]['protein_score']= float(row["RDS"])
                log.info(
                    'contig RDS: contig=%s, RDS=%f, score-sum=%f, #ORFs=%d',
                    contig_id, float(row["RDS"]), float(row["score_sum"]), int(row["ORFs"]))
                  
    scored_contigs = {}
    if(args.characterize):
        scored_contigs = contigs
    else:
        if(args.verbose):
            print(f'apply RDS sensitivity threshold (SNT={pc.RDS_SENSITIVITY_THRESHOLD:2.1f}) filter...')
            
        for (k, v) in contigs.items(): 
            if 'protein_score' not in v.keys():
                v['protein_score'] = 0
                scored_contigs[k] = v
            elif v['protein_score'] >= pc.RDS_SENSITIVITY_THRESHOLD:
                scored_contigs[k] = v
           
        no_excluded_contigs = len(contigs) - len(scored_contigs)
        log.info('RDS SNT filter: # discarded contigs=%d, # remaining contigs=%d' % (no_excluded_contigs, len(scored_contigs)))
        if(args.verbose):
            print(f'\texcluded {no_excluded_contigs} contigs by SNT filter')
            print('characterize contigs...')
            
    # the first Snakemake session starts for non-characsteristic part. 
    name = pf.get_base_name(cfg.genome_path)
    
    pf.fasta_into_chunk_contigs(scored_contigs, pc.DEFAULT_CONTIG_SIZE, output_path)
    pf.faa_into_chunk_contigs(pc.DEFAULT_CONTIG_SIZE, output_path)
      
    # second snakemake part starts for the characteristic part.
    result2 = snakemake.snakemake(
        snakefile="snakefiles/Snakefile2",
        workdir= output_path,
        config={"current_path": Path.cwd(), "min_circ_basepair_overlap":pc.MIN_CIRC_BASEPAIR_OVERLAP},
        cores=cfg.threads, scheduler="greedy", quiet=True)
    
    
    if result2:
        if(args.verbose):
            print("Workflow for characterization completed successfully!")
        log.info("Workflow for characterization completed.")
    else:
        if(args.verbose):
            print("Workflow for characterization failed.")
        log.debug("Workflow for characterization failed.")

    # parse the result from the characteristic part.
    hmm = ['amr_hits', 'conjugation_hits', 'mobilization_hits', 'replication_hits']                  

    for names in [ ['amr.tsv','amr_hits'], ['rrnas.tsv','rrnas'], ['orit.tsv','orit_hits'], ['cir.tsv','is_circular'], ['inc.tsv','inc_types'], 
            ['ref.tsv','plasmid_hits'], ['mob.tsv','mobilization_hits'], ['rep.tsv','replication_hits'], ['conj.tsv','conjugation_hits']]:
            
        if names[1] in hmm: 
            pf.extract_function_info_hmm(scored_contigs, names[0], names[1], output_path)
        elif names[1] == 'rrnas':
            pf.extract_function_info_rnnas(scored_contigs, names[0], names[1], output_path)
        elif names[1] == 'inc_types':
            pf.extract_function_info_inc(scored_contigs, names[0], names[1], output_path)
        elif names[1] == 'plasmid_hits':
            pf.extract_function_info_ref(scored_contigs, names[0], names[1], output_path)
        elif names[1] == 'orit_hits':
            pf.extract_function_info_orit(scored_contigs, names[0], names[1], output_path)
        elif names[1] == 'is_circular':
            pf.extract_function_info_cir(scored_contigs, names[0], names[1], output_path)
        else:
            pf.extract_function_info(scored_contigs, names[0], names[1], output_path)

    # lookup AMR genes with creating amr_genes dictionary and then match hmm-id.
    amr_genes = {}
    with cfg.db_path.joinpath('ncbifam-amr.tsv').open() as fh:
        for line in fh:
            cols = line.split('\t')
            amr_genes[cols[0]] = {
                'gene': cols[4],
                'product': cols[8]
            }
    
    for id, contig in scored_contigs.items():
        for hit in contig['amr_hits']:
            amr_gene = amr_genes[hit['hmm-id']]
            hit['gene'] = amr_gene['gene']
            hit['product'] = amr_gene['product']

    
    # filter contigs base on the mode: sensitivity, specificity, and accuracity.
    if(args.characterize):  
        filtered_contigs = scored_contigs
    elif(args.mode == 'sensitivity'):  
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if v['rrnas'] == 0}
    elif(args.mode == 'specificity'):
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if v['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD or v.get('protein_score', 0)} # fix 
    else:
        filtered_contigs = {k: v for (k, v) in scored_contigs.items() if pf.filter_contig(v)}


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

     # remove tmp dir
    shutil.rmtree(output_path.joinpath('tmp'), ignore_errors=True)
    log.debug('removed tmp dir: %s', cfg.tmp_path)

    return raw_contigs, contigs, filtered_contigs
