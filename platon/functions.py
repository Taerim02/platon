import logging
import subprocess as sp
import pyrodigal
import os, re
import pyhmmer
import pyfastx
import ast
import csv

import platon
import platon.config as cfg
import platon.constants as pc

log = logging.getLogger('functions')

def get_base_name(file_name):
    return os.path.splitext(os.path.basename(file_name))[0]

def write_sequence_to_file(contig, contig_split_position, header):
    file_path = cfg.tmp_path.joinpath(f"{contig['id']}-{header}.fasta")
    if header == 'a':
        sequence = contig['sequence'][:contig_split_position]
    elif header == 'b':
        sequence = contig['sequence'][contig_split_position:]
    with file_path.open(mode='w') as fh:
        fh.write('>{header}\n')
        fh.write(sequence +'\n ')
    return file_path, sequence

def run_command(cmd, env=False):
    env = cfg.env if env else None
    proc = sp.run(
        cmd,
        cwd=str(cfg.tmp_path),
        env=env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    return proc

def proc_error(var, var_error, proc, cmd, contig):
    if(proc.returncode != 0):
        log.warning(
            '%s failed! contig=%s, %s=%d',
            var, contig['id'], var_error, proc.returncode
        )
        log.debug(
            '%s: contig=%s, cmd=%s, stdout=\'%s\', stderr=\'%s\'',
            var, contig['id'], cmd, proc.stdout, proc.stderr
        )
        return


def contigs_into_chunks(contigs, contig_size, output_path):
    current_list = []
    current_size = 0
    i = 1
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    for id, contig in contigs.items():
        contig_name = contig['id']
        contig_sequence = contig['sequence']
        contig_length = contig['length']
        
        current_list.append({'contig_name': contig_name, 'contig_sequence': contig_sequence})
        current_size += contig_length
        
        if current_size >= contig_size:
            contig_file_path = os.path.join(f'{output_path}/tmp/chunk', f'{name}_{i}.fasta')
            with open(contig_file_path, 'a') as contig_file:
                for saved_contig in current_list:
                    contig_file.write(f">{saved_contig['contig_name']}\n")
                    contig_file.write(f"{saved_contig['contig_sequence']}\n")
            current_list = []
            current_size = 0
            i += 1 
    
    # Handle the last batch of contigs
    if len(current_list) > 0:
        i += 1
        contig_file_path = os.path.join(f'{output_path}/tmp/chunk', f'{name}_{i}.fasta')
        with open(contig_file_path, 'a') as contig_file:
            for saved_contig in current_list:
                contig_file.write(f">{saved_contig['contig_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")
    return
    
def fasta_into_chunk_contigs(contigs, contig_size, output_path):
    current_list = []
    current_size = 0
    i = 1
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    for id, contig in contigs.items():
        contig_name = contig['id']
        contig_sequence = contig['sequence']
        contig_length = contig['length']
        
        current_list.append({'contig_name': contig_name, 'contig_sequence': contig_sequence})
        current_size += contig_length
        
        if current_size >= contig_size:
            contig_file_path = os.path.join(f'{output_path}/tmp',  f'{name}_{i}_filtered.fasta')
            with open(contig_file_path, 'a') as contig_file:
                for saved_contig in current_list:
                    contig_file.write(f">{saved_contig['contig_name']}\n")
                    contig_file.write(f"{saved_contig['contig_sequence']}\n")
            current_list = []
            current_size = 0
            i += 1 
    
    # Handle the last batch of contigs
    if len(current_list) > 0:
        i += 1
        contig_file_path = os.path.join(f'{output_path}/tmp',  f'{name}_{i}_filtered.fasta')
        with open(contig_file_path, 'a') as contig_file:
            for saved_contig in current_list:
                contig_file.write(f">{saved_contig['contig_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")
    return


def faa_into_chunk_contigs(contig_size, output_path):
    contigs_to_save = []
    current_size = 0
    i = 0
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    protein_files = os.listdir(output_path.joinpath(f'{output_path}/tmp/protein'))
    for file in protein_files:
      if file.endswith(".faa"):
          for record in pyfastx.Fasta(os.path.join(output_path.joinpath(f'{output_path}/tmp/protein'), file)):
              orf_name = str(record.name).split()[0]
              contig_sequence = str(record.seq)
              contig_length = len(record.seq)
              
              if current_size + contig_length <= contig_size:
                  contigs_to_save.append({'orf_name': orf_name, 'contig_sequence': contig_sequence})
                  current_size += contig_length
              else:
                  i += 1
                  contig_file_path = os.path.join(f'{output_path}/tmp', f'{name}_{i}_filtered.faa')
                  with open(contig_file_path, 'w') as contig_file:
                      for saved_contig in contigs_to_save:
                          contig_file.write(f">{saved_contig['orf_name']}\n")
                          contig_file.write(f"{saved_contig['contig_sequence']}\n")
                  contigs_to_save = []
                  current_size = 0
    
    # Handle the last batch of contigs
    if contigs_to_save:
        i += 1
        contig_file_path = os.path.join(f'{output_path}/tmp', f'{name}_{i}_filtered.faa')
        with open(contig_file_path, 'w') as contig_file:
            for saved_contig in contigs_to_save:
                contig_file.write(f">{saved_contig['orf_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")
    return
                
def test_circularity(contig):
    """Test if this contig can be circularized."""
    contig_split_position = int(contig['length'] / 2)
    contig_fragment_a_path, contig_fragment_a_seq = write_sequence_to_file(contig, contig_split_position, 'a')
    contig_fragment_b_path, contig_fragment_b_seq = write_sequence_to_file(contig, contig_split_position, 'b')
    log.debug(
        'circularity: contig=%s, len=%d, seq-a-len=%d, seq-b-len=%d',
        contig['id'], contig['length'], len(contig_fragment_a_seq), len(contig_fragment_b_seq)
    )
    cmd = [
        'nucmer',
        '-f',  # only forward strand
        '-l', '40',  # increase min match length to 40 bp
        '--threads=1',
        '-p', contig['id'],
        str(contig_fragment_b_path),
        str(contig_fragment_a_path)
    ]
    proc = run_command(cmd, env=False)
    proc_error('circularity', 'nucmer-error-code', proc, cmd, contig)
    has_match = False
    with cfg.tmp_path.joinpath(f"{contig['id']}.delta").open() as fh:
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
                            and alignment_a > pc.MIN_CIRC_BASEPAIR_OVERLAP
                            and (mismatches / alignment_a) < 0.05
                            and end_b == len(contig_fragment_b_seq)
                            and start_a == 1):
                        contig['is_circular'] = True
                        link = {
                            'length': alignment_a,
                            'mismatches': mismatches,
                            'prime5Start': 1,
                            'prime5End': alignment_a,
                            'prime3Start': contig['length'] - alignment_b + 1,
                            'prime3End': contig['length'],
                        }
                        contig['circular_link'] = link
                        log.info(
                            'circularity: link! id=%s, length=%d, # mismatches=%d, linking-region=[1-%d]...[%d-%d]',
                            contig['id'], link['length'], link['mismatches'], link['prime5End'], link['prime3Start'], link['prime3End']
                        )
                        break
    log.info('circularity: contig=%s, is-circ=%s', contig['id'], contig['is_circular'])
    return


def search_inc_type(contig):    ## why they do not record for the coverage and identity?
    """Search for incompatibility motifs."""
    contig_path = cfg.tmp_path.joinpath(f"{contig['id']}.fasta")
    tmp_output_path = cfg.tmp_path.joinpath(f"{contig['id']}.inc.blast.out")
    cmd = [
        'blastn',
        '-query', str(cfg.db_path.joinpath('inc-types.fasta')),
        '-subject', str(contig_path),
        '-num_threads', '1',
        '-perc_identity', '90',
        '-culling_limit', '1',
        '-outfmt', '6 qseqid sstart send sstrand pident qcovs bitscore',
        '-out', str(tmp_output_path)
    ]
    proc = run_command(cmd, env=True)
    proc_error('inc-type', 'blastn-error-code', proc, cmd, contig)
    hits_per_pos = {}
    with tmp_output_path.open() as fh:
        for line in fh:
            cols = line.rstrip().split('\t')
            #log.info('inc:', len(cols))
            hit = {
                'type': cols[0],
                'start': int(cols[1]),
                'end': int(cols[2]),
                'strand': '+' if cols[3] == 'plus' else '-',
                'identity': float(cols[4]) / 100,
                'coverage': float(cols[5]) / 100,
                'bitscore': int(cols[6])
            }
            if(hit['coverage'] >= 0.6):
                hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                if(hit_pos in hits_per_pos):
                    former_hit = hits_per_pos[hit_pos]
                    if(hit['bitscore'] > former_hit['bitscore']):
                        hits_per_pos[hit_pos] = hit
                        log.info(
                            'inc-type: hit! contig=%s, inc id=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
                else:
                    hits_per_pos[hit_pos] = hit
                    log.info(
                        'inc-type: hit! contig=%s, inc id=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    contig['inc_types'] = list(hits_per_pos.values())
    log.info('inc-type: contig=%s, # inc-types=%s', contig['id'], len(contig['inc_types']))
    return


def search_rrnas(contig):
    """Search for ribosomal RNA sequences."""
    contig_path = cfg.tmp_path.joinpath(f"{contig['id']}.fasta")
    tmp_output_path = cfg.tmp_path.joinpath(f"{contig['id']}.rrna.cmscan.tsv")
    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '--cpu', '1',
        '--tblout', str(tmp_output_path),
        str(cfg.db_path.joinpath('rRNA')),
        str(contig_path)
    ]
    proc = run_command(cmd, env=False)
    proc_error('rRNAs', 'cmscan-error-code', proc, cmd, contig)
    with tmp_output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                cols = line.rstrip().split()
                #log.info('rnnas:', cols)
                hit = {
                    'type': cols[0],
                    'start': int(cols[7]),
                    'end': int(cols[8]),
                    'strand': cols[9],
                    'bitscore': float(cols[14]),
                    'evalue': float(cols[15])
                }
                contig['rrnas'].append(hit)
                log.info(
                    'rRNAs: hit! contig=%s, rRNA id=%s, start=%d, end=%d, strand=%s',
                    contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                )
    log.info('rRNAs: contig=%s, # rRNAs=%s', contig['id'], len(contig['rrnas']))
    return



def search_reference_plasmids(contig):
    """Search for reference plasmid hits."""

    # Reduce blastn word size to overcome segmentation faults due to too many HSPs.
    blast_word_size = int(contig['length'] / 1000)
    if(blast_word_size < 11):
        blast_word_size = 11

    contig_path = cfg.tmp_path.joinpath(f"{contig['id']}.fasta")
    tmp_output_path = cfg.tmp_path.joinpath(f"{contig['id']}.refplas.blast.out")

    cmd = [
        'blastn',
        '-query', str(contig_path),
        '-db', str(cfg.db_path.joinpath('refseq-plasmids')),
        '-num_threads', '1',
        '-culling_limit', '1',
        '-perc_identity', '80',
        '-word_size', str(blast_word_size),
        '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
        '-out', str(tmp_output_path)
    ]
    proc = run_command(cmd, env=True)
    proc_error('ref plasmids', 'blastn-error-code', proc, cmd, contig)
    with tmp_output_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split('\t')
            hit = {
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'plasmid_start': int(cols[3]),
                'plasmid_end': int(cols[4]),
                'plasmid': {
                    'id': cols[0],
                    'length': int(cols[5])
                },
                'coverage': float(cols[6]) / contig['length'],
                'identity': float(cols[7]) / float(cols[6])
            }
            if(hit['coverage'] >= 0.8):
                contig['plasmid_hits'].append(hit)
                log.info(
                    'ref plasmids: hit! contig=%s, ref id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    contig['id'], hit['plasmid']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('ref plasmids: contig=%s, # ref plasmids=%s', contig['id'], len(contig['plasmid_hits']))
    return


def search_orit_sequences(contig):
    """Search for oriT sequence hits."""

    contig_path = cfg.tmp_path.joinpath(f"{contig['id']}.fasta")
    tmp_output_path = cfg.tmp_path.joinpath(f"{contig['id']}.orit.blast.out")
    cmd = [
        'blastn',
        '-query', str(contig_path),
        '-db', str(cfg.db_path.joinpath('orit')),
        '-num_threads', '1',
        '-culling_limit', '1',
        '-perc_identity', '90',
        '-evalue', '1E-5',
        '-outfmt', '6 sseqid qstart qend sstart send slen length nident',
        '-out', str(tmp_output_path)
    ]
    proc = run_command(cmd, env=True)
    proc_error('oriT', 'blastn-error-code', proc, cmd, contig)
    with tmp_output_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            cols = line.split('\t')
            hit = {
                'contig_start': int(cols[1]),
                'contig_end': int(cols[2]),
                'orit_start': int(cols[3]),
                'orit_end': int(cols[4]),
                'orit': {
                    'id': cols[0],
                    'length': int(cols[5])
                },
                'coverage': float(cols[6]) / int(cols[5]),
                'identity': float(cols[7]) / float(cols[6])
            }
            if(hit['coverage'] >= 0.9 and hit['identity'] >= 0.9):
                contig['orit_hits'].append(hit)
                log.info(
                    'oriT: hit! contig=%s, oriT id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    contig['id'], hit['orit']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('oriT: contig=%s, # oriT=%s', contig['id'], len(contig['orit_hits']))
    return
        
        
def filter_contig(contig):
    
    """Apply heuristic filters based on contig information."""
    # include all circular contigs
    if(contig['is_circular']):
        log.debug('filter: is circ! contig=%s', contig['id'])
        return True

    # include all contigs with Inc type signatures
    if(len(contig['inc_types']) > 0):
        log.debug('filter: has inc types! contig=%s', contig['id'])
        return True

    # include all contigs containing replication genes
    if(len(contig['replication_hits']) > 0):
        log.debug('filter: has rep hits! contig=%s', contig['id'])
        return True

    # include all contigs containing mobilization genes
    if(len(contig['mobilization_hits']) > 0):
        log.debug('filter: has mob hits! contig=%s', contig['id'])
        return True

    # include all contigs containing oriT sequences
    if(len(contig['orit_hits']) > 0):
        log.debug('filter: has oriT hits! contig=%s', contig['id'])
        return True

    # include all contigs with high confidence protein scores
    if('protein_score' in contig.keys() and contig['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD):
        log.debug('filter: RDS > SPT! contig=%s', contig['id'])
        return True

    # include all contigs with mediocre protein scores but additional blast hit evidence without rRNAs
    if('protein_score' in contig.keys() and contig['protein_score'] >= pc.RDS_CONSERVATIVE_THRESHOLD 
        and len(contig['plasmid_hits']) > 0
        and len(contig['rrnas']) == 0):
        log.debug('filter: RDS > CT & plasmid hits & no rRNAs! contig=%s', contig['id'])
        return True
            
    return False
    

def train_gene_prediction(contigs): 
    """ Train gene prodiction for the genomic mode"""
    log.info('create prodigal training info object: meta=%s, not pyrodigal_metamode')
    training_info = None
    orf_finder = pyrodigal.GeneFinder(meta= False, closed=True)
    seqs = [contig['sequence'] for id, contig in contigs.items()]
    combined_sequence = ''.join(seqs)
    bytes_combined_sequence = combined_sequence.encode('utf-8')
    training_info = orf_finder.train(bytes_combined_sequence) ##  list comprehsion
    return training_info
    

def predict_orfs_py(contigs, record, proteins_path, pyrodigal_metamode, training_info = None):
    orf_finder = pyrodigal.GeneFinder(training_info, meta=pyrodigal_metamode, closed=True, mask=True)
    genes = orf_finder.find_genes(str(record.seq))
    with open(str(proteins_path), "a") as dst:
        genes.write_translations(dst, sequence_id=record.name)
    for i, gene in enumerate(genes):
            orf_id = i+1
            strand_match = re.search(r"strand=([-+])", str(gene))
            orf = {
                'start': gene.begin,
                'end': gene.end,
                'strand': strand_match.group(1),
                'id': orf_id
            }
            contig = contigs.get(record.name, None)    
            if(contig is not None):
                contig['orfs'][orf_id] = orf
                log.info(
                    'ORFs: found! contig=%s, start=%d, end=%d, strand=%s',
                    contig['id'], orf['start'], orf['end'], orf['strand']
                )
    return 
    

def construct_hit_dict(hits, hit, orf, inclued_hmm_id = False):
    if inclued_hmm_id:
      hmm_id = hit.best_domain.alignment.hmm_accession.decode()
    hit = {
        'type': hits.query_name.decode(),
        'start': int(orf['start']),
        'end': int(orf['end']),
        'strand': orf['strand'],
        'bitscore': hit.score,
        'evalue': hit.evalue
    }
    if inclued_hmm_id:
        hit['hmm-id'] = hmm_id
    return hit
    

def search_amr_genes_py(contigs, filteredProteinsPath):
    """Search for AMR genes."""
    hits_set = set()
    hmm_profile_path = str(cfg.db_path.joinpath('ncbifam-amr'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, bit_cutoffs="trusted"):
            for hit in hits:
                if hit.included:
                    contig_id, orf_id = hit.name.decode().rsplit('_', 1)
                    if f'{contig_id}_{orf_id}' not in hits_set:
                        contig = contigs[contig_id]
                        orf = contig['orfs'][int(orf_id)]
                        hit = construct_hit_dict(hits, hit, orf, inclued_hmm_id = True)
                        hits_set.add(f'{contig_id}_{orf_id}')
                        contig['amr_hits'].append(hit)
                        log.info(
                            'AMRs: hit! contig=%s, AMR id=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
    for id, contig in contigs: 
        log.info('AMRs: contig=%s, # AMRs=%s', contig['id'], len(contig['amr_hits']))
    return


def search_replication_genes_py(contigs, filteredProteinsPath):
    """Search for replication genes."""
    hits_set = set()
    hmm_profile_path = str(cfg.db_path.joinpath('replication'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-100):
            for hit in hits:
                if hit.included:
                    contig_id, orf_id = hit.name.decode().rsplit('_', 1)
                    if f'{contig_id}_{orf_id}' not in hits_set:
                        contig = contigs[contig_id]
                        orf = contig['orfs'][int(orf_id)]
                        hit = construct_hit_dict(hits, hit, orf)
                        hits_set.add(f'{contig_id}_{orf_id}')
                        contig['replication_hits'].append(hit)
                        log.info(
                            'rep genes: hit! contig=%s, rep id=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
    for id, contig in contigs:                
        log.info('rep genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['replication_hits']))
    return

def search_mobilization_genes_py(contigs, filteredProteinsPath):
    """Search for mobilization genes."""
    hits_set = set()
    hmm_profile_path = str(cfg.db_path.joinpath('mobilization'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-10):
            for hit in hits:
                if hit.included:
                    contig_id, orf_id = hit.name.decode().rsplit('_', 1)
                    if f'{contig_id}_{orf_id}' not in hits_set:
                        contig = contigs[contig_id]
                        orf = contig['orfs'][int(orf_id)]
                        hit = construct_hit_dict(hits, hit, orf)
                        hits_set.add(f'{contig_id}_{orf_id}')
                        contig['mobilization_hits'].append(hit)
                        log.info(
                            ' mob genes: hit! contig=%s, mob id=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
    for id, contig in contigs:                
        log.info('mob genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['mobilization_hits']))
    return

def search_conjugation_genes_py(contigs, filteredProteinsPath):
    """Search for conjugation genes."""
    hits_set = set()
    hmm_profile_path = str(cfg.db_path.joinpath('conjugation'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-100):
            for hit in hits:
                if hit.included:
                    contig_id, orf_id = hit.name.decode().rsplit('_', 1)
                    if f'{contig_id}_{orf_id}' not in hits_set:
                        contig = contigs[contig_id]
                        orf = contig['orfs'][int(orf_id)]
                        hit = construct_hit_dict(hits, hit, orf)
                        hits_set.add(f'{contig_id}_{orf_id}')
                        contig['conjugation_hits'].append(hit)
                        log.info(
                            'conj genes: hit! contig=%s, conj id=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
    for id, contig in contigs:
        log.info('conj genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['conjugation_hits']))
    return
    
def merge_dicts(dict1, dict2):
    return {**dict1, **dict2}
    
def extract_function_info_rnnas(contigs:dict, tsv_file:str, index:str, output_path): 
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row['contig'] in contigs:
                hit = {
                    'type': row['type'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'],
                    'bitscore': float(row['bitscore']),
                    'evalue': float(row['evalue'])
                }
                contigs[row['contig']][index].append(hit)
                log.debug(
                    'rRNAs: hit! contig=%s, rRNA id=%s, start=%d, end=%d, strand=%s',
                    row["contig"], hit['type'],hit['start'], hit['end'], hit['strand']
                )
    log.info('%s: found %d number contigs!', index, len([k for k, v in contigs.items() if len(v[index]) > 0]))
    return 
        

def extract_function_info_inc(contigs:dict, tsv_file:str, index:str, output_path): 
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t") 
        hits_dict = {}
        for row in reader:
            if row['contig'] in contigs:
                hit = {
                    'type': row['type'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'],
                    'identity': float(row['identity']),
                    'coverage': float(row['coverage']),
                    'bitscore': float(row['bitscore'])
                }
                if(hit['coverage'] >= 0.6):
                    hit_pos = hit['end'] if hit['strand'] == '+' else hit['start']
                    hit['hit_pos'] = hit_pos
                    if row['contig'] in hits_dict:
                        hits_dict[row['contig']][hit_pos] = hit
                    else:
                        hits_dict[row['contig']] = {hit_pos:hit}
                    log.debug(
                        'inc hit: hit! contig=%s, inc id=%s, start=%d, end=%d, strand=%s',
                        row["contig"], hit['type'],hit['start'], hit['end'], hit['strand']
                    )
    for key, value in hits_dict.items():
        contigs[key]['inc_types'] = list(value.values())
    log.info('%s: found %d number contigs!', index, len(list(hits_dict.keys())))
    return 
    
    
def extract_function_info_ref(contigs:dict, tsv_file:str, index:str, output_path): 
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t") 
        for row in reader:
            if row['contig'] in contigs:
                hit = {
                'contig_start': int(row['contig_start']),
                'contig_end': int(row['contig_end']),
                'plasmid_start': int(row['plasmid_start']),
                'plasmid_end': int(row['plasmid_end']),
                'plasmid': {
                    'id': row['plasmid_id'],
                    'length': int(row['plasmid_length'])
                },
                'coverage': float(row['coverage']),
                'identity': float(row['identity'])
                }
                contigs[row['contig']]['plasmid_hits'].append(hit)
                log.debug(
                    'ref plasmids: hit! contig=%s, ref id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    row['contig'], hit['plasmid']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('%s: found %d number contigs!', index, len([k for k, v in contigs.items() if len(v[index]) > 0])) 
    return
    
def extract_function_info_orit(contigs:dict, tsv_file:str, index:str, output_path): 
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t") 
        for row in reader:
            if row['contig'] in contigs:
                hit = {
                'contig_start': int(row['contig_start']),
                'contig_end': int(row['contig_end']),
                'orit_start': int(row['orit_start']),
                'orit_end': int(row['orit_end']),
                'orit': {
                    'id': row['orit_id'],
                    'length': int(row['orit_length'])
                },
                'coverage': float(row['coverage']),
                'identity': float(row['identity'])
                }
                contigs[row['contig']]['orit_hits'].append(hit)
                log.debug(
                    'oriT plasmids: hit! contig=%s, oriT id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    row['contig'], hit['orit']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('%s: found %d number contigs!', index, len([k for k, v in contigs.items() if len(v[index]) > 0]))  
    return
    

def extract_function_info_cir(contigs:dict, tsv_file:str, index:str, output_path): 
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t") 
        for row in reader:
            if row['contig'] in contigs:
                link = {
                    'length': row['length'],
                    'mismatches': row['mismatches'],
                    'prime5Start': 1,
                    'prime5End': row['prime5End'],
                    'prime3Start': row['prime3Start'],
                    'prime3End': row['prime3End'],
                        }
                contigs[row['contig']]['is_circular'] = True
                contigs[row['contig']]['circular_link'] = link
                log.debug('circularity: contig=%s, len=%d, seq-a-len=%d, seq-b-len=%d',
                                row["contig"], contigs[row['contig']]['length'], row['seq-a-len'], row['seq-b-len'])
                log.debug('circularity: link! id=%s, length=%d, # mismatches=%d, linking-region=[1-%d]...[%d-%d]',
                            row["contig"], link['length'], link['mismatches'], link['prime5End'], link['prime3Start'], link['prime3End'])
    log.info('%s: found %d number contigs!', index, len([k for k, v in contigs.items() if v[index]]))  
    return
                          
def extract_function_info_hmm(contigs:dict, tsv_file:str, index:str, output_path:str):         
    hits_set = set()
    with open(os.path.join(output_path.joinpath('tmp/function'), tsv_file),"r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row['contig'] in contigs:
                if(f'{row["contig"]}_{row["orf_id"]}' not in hits_set):
                    orf = contigs[row["contig"]]['orfs'][int(row["orf_id"])]
                    hit = {
                        'type': row['type'],
                        'start': int(orf['start']),
                        'end': int(orf['end']),
                        'strand': orf['strand'],
                        'bitscore': float(row['bitscore']),
                        'evalue': float(row['evalue'])
                    }
                    if index == 'amr_hits':
                        hit['hmm-id'] = row['hmm-id']
                    hits_set.add(f'{row["contig"]}_{row["orf_id"]}')
                    contigs[row['contig']][index].append(hit)
                    index_dict = {'amr_hits':['AMRs', 'AMR'], 'conjugation_hits': ['conj genes', 'conj'], 'mobilization_hits': ['mob genes', 'mob'], 'replication_hits': ['rep genes', 'rep']}
                    log.debug(
                        '%s: hit! contig=%s, %s id=%s, start=%d, end=%d, strand=%s',
                        index_dict[index][0], row["contig"], index_dict[index][1], hit['type'],hit['start'], hit['end'], hit['strand']
                    )
    log.info('%s: found %d number contigs!', index, len([k for k, v in contigs.items() if len(v[index]) > 0]))                
    return
