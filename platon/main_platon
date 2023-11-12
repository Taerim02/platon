import db as db
import config as cfg
import constants as pc
import functions as pf
import utils as pu
import meta_platon
import single_platon
from pathlib import Path
import cProfile
import functools as ft
import json
import logging
import os, re, sys
import shutil
import subprocess as sp
import __init__ as init

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

    if cfg.module == 'metagenomic':
        meta_platon.main()
    else:
        single_platon.main()


if __name__ == '__main__':
    main()
    #profiler = cProfile.Profile()
    #profiler.runctx("main()", globals(), locals())
