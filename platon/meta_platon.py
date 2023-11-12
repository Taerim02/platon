import functools as ft
import json
import logging
import os
import sys
from pathlib import Path
import pyfastx
import cProfile

import __init__ as init
import db as db
import config as cfg
import constants as pc
import functions as pf
import utils as pu
from snakemake import snakemake
import pandas as pd
import re
import ast


def main():
    # parse arguments

    args = pu.parse_arguments()
  
    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.genome).stem
    
    try:
        output_path = Path(args.output) if args.output else Path.cwd()
        if(not output_path.exists()):
            output_path.mkdir(parents=True, exist_ok=True)
        elif(not os.access(str(output_path), os.X_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not accessible!')
        elif(not os.access(str(output_path), os.W_OK)):
            sys.exit(f'ERROR: output path ({output_path}) not writable!')
        output_path = output_path.resolve()
        cfg.output_path = output_path
        
    except:
        sys.exit(f'ERROR: could not resolve or create output directory ({args.output})!')
    logging.basicConfig(
        filename=str(output_path.joinpath(f'{cfg.prefix}.log')),
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO
    )
    log = logging.getLogger('MAIN')
    log.info('version %s', init.__version__)
    log.info('command line: %s', ' '.join(sys.argv))

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    db.check()  # test database
    pu.test_dependencies()  # test dependencies
    
    if(cfg.verbose):
        print(f'Platon v{init.__version__}')
        print('Options and arguments:')
        print(f'\tmain input: {cfg.genome_path}')
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}")
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tmode: {cfg.mode}')
        print(f'\tcharacterize: {cfg.characterize}')
        print(f'\ttmp path: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')

    # parse draft genome
    if(args.verbose):
        print('parse draft genome...')
    contigs = {}
    raw_contigs = []
    try:
        for record in pyfastx.Fasta(str(cfg.genome_path)):
            length = len(record.seq)
            contig = {
                'id': record.name,
                'length': length,
                'sequence': str(record.seq),  # is there any other way to keep this seq?
                'orfs': {},
                'is_circular': False,
                'inc_types': [],
                'amr_hits': [],
                'mobilization_hits': [],
                'orit_hits': [],
                'replication_hits': [],
                'conjugation_hits': [],
                'rrnas': [],
                'plasmid_hits': []
            }
            raw_contigs.append(contig)
                       
            match_spades = re.fullmatch(pc.SPADES_CONTIG_PATTERN, record.name)
            match_unicycler = re.fullmatch(pc.UNICYCLER_CONTIG_PATTERN, record.description)
            if(match_spades is not None):
                contig['coverage'] = float(match_spades.group(1))
            elif(match_unicycler is not None):
                contig['coverage'] = float(match_unicycler.group(1))
                if(match_unicycler.group(2) is not None):
                    contig['is_circular'] = True                  # can be circularized
                    contig['type'] = True
            else:
                contig['coverage'] = 0
            
            # only include contigs with reasonable lengths except of platon runs in characterization mode
            
            if(args.characterize):
                contigs[record.name] = contig
            else:                                            
                if(length < pc.MIN_CONTIG_LENGTH):
                    log.info('exclude contig: too short: id=%s, length=%d', record.name, length)
                    if (args.verbose):
                        print(f'\texclude contig \'{record.name}\', too short ({length})')
                elif(length >= pc.MAX_CONTIG_LENGTH):
                    log.info('exclude contig: too long: id=%s, length=%d', record.name, length)
                    if (args.verbose):
                        print(f'\texclude contig \'{record.name}\', too long ({length})')
                else:
                    contigs[record.name] = contig

    except: 
        log.error('something went wrong!', exc_info=True)
        sys.exit('ERROR: something went wrong!')
        
    for name in ['rrnas', 'orit', 'inc', 'ref', 'orf', 'function']:
        if not os.path.exists(f'tmp/{name}'):
            os.mkdir(f'tmp/{name}')


    i_list = pf.contigs_into_chunks(contigs, pc.DEFAULT_CONTIG_SIZE) 

    result1 = snakemake(
        snakefile="Snakefile_meta_1.py",
        workdir=cfg.output_path,
        cores=10, 
    )

    pf.delete_contig_files(i_list)

    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    
    full_filtered_contig_path = os.path.join('tmp', f"{name}_filtered.fasta")
    pf.fasta_into_chunk(full_filtered_contig_path, pc.DEFAULT_CONTIG_SIZE)
    os.remove(full_filtered_contig_path)

    full_filtered_proteins_path = os.path.join('tmp', f"{name}_filtered.faa")
    pf.faa_into_chunk(full_filtered_proteins_path, pc.DEFAULT_CONTIG_SIZE)
    os.remove(full_filtered_proteins_path)

    result2 = snakemake(
        snakefile="Snakefile_meta_2_blast.py",
        workdir=cfg.output_path,
        cores=20,  scheduler="greedy" # keepgoing=True
    )
    
    
    orf_files = os.listdir('tmp/orf')
    for file in orf_files:
        with open(os.path.join('tmp/orf', file),"r") as fh:
            for line in fh:
                line = json.loads(line) 
                contig_id = line['contig']
                orf_id = line['id']
                contigs[contig_id]['orfs'][orf_id]= line
    
    
    protein_score_pattern = r'contig RDS: contig=([A-Z0-9]+), RDS=(-?[0-9.]+)'
    protein_score_files = os.listdir('tmp/protein_score')
    for file in protein_score_files:
        with open(os.path.join('tmp/protein_score', file),"r") as fh:
            for line in fh:
                match = re.match(protein_score_pattern, line)
                if match:
                    #if contig_id in contig_ids:  
                    contig_id = match.group(1)
                    contigs[contig_id]['protein_score']= float(match.group(2))
                    
    for name_dict in [{'amr.txt':'amr_hits'}, {'rrnas.txt':'rrnas'}, {'orit.txt':'orit_hits'}, {'cir.txt':'is_circular'}, {'inc.txt':'inc_types'}, 
            {'ref.txt':'plasmid_hits'}, {'mob.txt': 'mobilization_hits'}, {'rep.txt':'replication_hits'}, {'conj.txt':'conjugation_hits'}]:
            pf.extract_function_info(contigs, list(name_dict.keys())[0], name_dict[list(name_dict.keys())[0]])
                      
    """                
    function_pattern = r"id: ([A-Z0-9]+) \{(.+)\}"

    with open(os.path.join('tmp/function', 'amr.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:  
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['amr_hits'].append(dict_value)


    with open(os.path.join('tmp/function', 'rrnas.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['rrnas'].append(match.group(2))

    with open(os.path.join('tmp/function', 'orit.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:  
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['orit_hits'].append(dict_value)
    

    with open(os.path.join('tmp/function', 'cir.txt'),"r") as fh: 
        for line in fh:
            match = re.match(function_pattern, line)
            if match:
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['is_circular'].append(match.group(2))
                    

    with open(os.path.join('tmp/function', 'inc.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:  
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['inc_types'].append(dict_value)
                    

    with open(os.path.join('tmp/function', 'ref.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:  
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['plasmid_hits'].append(dict_value)
                    
    with open(os.path.join('tmp/function', 'mob.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['mobilization_hits'].append(dict_value)
            

    with open(os.path.join('tmp/function', 'rep.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match:  
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                contigs[contig_id]['replication_hits'].append(dict_value)
                    
    with open(os.path.join('tmp/function', 'conj.txt'),"r") as fh:
        for line in fh:
            match = re.match(function_pattern, line)
            if match: 
                contig_id = match.group(1)
                dict_string = "{" + match.group(2) + "}"
                dict_value = ast.literal_eval(dict_string)
                print(dict_string)
                contigs[contig_id]['conjugation_hits'].append(dict_value)
    """
                    
    filtered_contigs = None
    if(args.characterize):  # skip protein score based filtering
        filtered_contigs = contigs
    elif(args.mode == 'sensitivity'):  # skip protein score based filtering but apply rRNA filter
        filtered_contigs = {k: v for (k, v) in contigs.items() if v['rrnas'] == 0}
    elif(args.mode == 'specificity'):
        filtered_contigs = {k: v for (k, v) in contigs.items() if v['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD}
    else:
        filtered_contigs = {k: v for (k, v) in contigs.items() if pf.filter_contig_meta(v)}
        
        
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
                  

    # write comprehensive results to JSON file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.json')
    log.debug('output: json=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        indent = '\t' if args.verbose else None
        separators = (', ', ': ') if args.verbose else (',', ':')
        json.dump(filtered_contigs, fh, indent=indent, separators=separators)

    # write chromosome contigs to fasta file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.chromosome.fasta')
    log.debug('output: chromosomes=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] not in filtered_contigs):
                fh.write(f">{contig['id']}\n{contig['sequence']}\n")

    # write plasmid contigs to fasta file
    tmp_output_path = output_path.joinpath(f'{cfg.prefix}.plasmid.fasta')
    log.debug('output: plasmids=%s', tmp_output_path)
    with tmp_output_path.open(mode='w') as fh:
        for contig in raw_contigs:
            if(contig['id'] in filtered_contigs):
                fh.write(f">{contig['id']}\n{contig['sequence']}\n")
            
        
    

if __name__ == '__main__':
    main()
    
