import json
import logging
import os
import re
import sys
from pathlib import Path

import pyfastx

import platon
import platon.config as cfg
import platon.constants as pc
import platon.db as db
import platon.utils as pu
import platon.__init__ as init


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
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}")
        print(f'\toutput: {cfg.output_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tmode: {cfg.mode}')
        print(f'\tcharacterize: {cfg.characterize}')
        print(f'\ttmp path: {cfg.tmp_path}')
        print(f'\t# threads: {cfg.threads}')
    
    # set lists for the JSON file.    
    raw_contigs = []
    filtered_contigs = None
    contigs = {}

    # parse draft genome
    if(args.verbose):
        print('parse draft genome...')
    try:
        for record in pyfastx.Fasta(str(cfg.genome_path)):
            length = len(record.seq)
            contig = {
                'id': record.name,
                'length': length,
                'sequence': str(record.seq),
                'orfs': {},
                'is_circular': False,
                'inc_types': [],
                'amr_hits': [],
                'mobilization_hits': [],
                'orit_hits': [],
                'replication_hits': [],
                'conjugation_hits': [],
                'rrnas': [],
                'plasmid_hits': [],
            }
            raw_contigs.append(contig)

            # read coverage from contig names if they were assembled with SPAdes
            match_spades = re.fullmatch(pc.SPADES_CONTIG_PATTERN, record.name)
            match_unicycler = re.fullmatch(pc.UNICYCLER_CONTIG_PATTERN, record.description)
            if(match_spades is not None):
                contig['coverage'] = float(match_spades.group(1))
            elif(match_unicycler is not None):
                contig['coverage'] = float(match_unicycler.group(1))
                if(match_unicycler.group(2) is not None):
                    contig['is_circular'] = True                  # can be circularized
            else:
                contig['coverage'] = 0

            # only include contigs with reasonable lengths except of
            # platon runs in characterization mode
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

    if(args.verbose):
        print(f'\tparsed {len(raw_contigs)} raw contigs')
        print(f'\texcluded {len(raw_contigs) - len(contigs)} contigs by size filter')
        print(f'\tanalyze {len(contigs)} contigs')
    log.info(
        'length contig filter: # input=%d, # discarded=%d, # remaining=%d',
        len(raw_contigs), (len(raw_contigs) - len(contigs)), len(contigs)
    )

    if(len(raw_contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')

    if(len(contigs) == 0):
        print(pc.HEADER)
        log.warning('no potential plasmid contigs!')
        if(args.verbose):
            print('No potential plasmid contigs found. Please, check contig lengths. Maybe you passed a finished or pseudo genome?')
        sys.exit(0)
    
    if cfg.metagenome:
        import platon.meta_platon as meta_platon
        raw_contigs, contigs, filtered_contigs = meta_platon.main(raw_contigs, contigs, args, log, output_path)
    else:
        import platon.single_platon as single_platon
        raw_contigs, contigs, filtered_contigs = single_platon.main(raw_contigs, contigs, args, log, output_path)
    
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
