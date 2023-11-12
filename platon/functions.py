import logging
import subprocess as sp

import config as cfg
import constants as pc
import pyrodigal
import re
import pyhmmer
import pandas as pd

import os 
import pyfastx

log = logging.getLogger('functions')

def get_base_name(file_name):
    return os.path.splitext(os.path.basename(file_name))[0]

def check_plasmid_true(contig, var_name):
    if contig[var_name] == True and not contig['type']: # avoid overwriting multiple times    
        contig['type'] = "plasmid"    
    return

def check_plasmid_length(contig, var_name):
    if len(contig[var_name]) > 0 and not contig['type']:
        contig['type'] = "plasmid"
    return

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

def delete_empty_file(file_path):
    if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
        os.remove(file_path)

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


def contigs_into_chunks(contigs, contig_size):
    contigs_to_save = []
    current_size = 0
    i_list = []
    i = 0
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    
    for contig in contigs.values():
        contig_name = contig['id']
        contig_sequence = contig['sequence']
        contig_length = contig['length']
        
        if current_size + contig_length <= contig_size:
            contigs_to_save.append({'contig_name': contig_name, 'contig_sequence': contig_sequence})
            current_size += contig_length
        else:
            i += 1
            i_list.append(i)
            contig_file_path = os.path.join(cfg.output_path, f'{name}_{i}.fasta')
            with open(contig_file_path, 'w') as contig_file:
                for saved_contig in contigs_to_save:
                    contig_file.write(f">{saved_contig['contig_name']}\n")
                    contig_file.write(f"{saved_contig['contig_sequence']}\n")
            contigs_to_save = []
            current_size = 0
    
    # Handle the last batch of contigs
    if contigs_to_save:
        i += 1
        contig_file_path = os.path.join(cfg.output_path, f'{name}_{i}.fasta')
        with open(contig_file_path, 'w') as contig_file:
            for saved_contig in contigs_to_save:
                contig_file.write(f">{saved_contig['contig_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")
    return i_list

def delete_contig_files(i_list):
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    for i in i_list:
        os.remove(os.path.join(cfg.output_path, f'{name}_{i}.fasta'))

def fasta_into_chunk(fasta_file, contig_size):
    contigs_to_save = []
    current_size = 0
    i = 0
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    for record in pyfastx.Fasta(str(fasta_file)):
        contig_name = record.name
        contig_sequence = str(record.seq)
        contig_length = len(record.seq)
        
        if current_size + contig_length <= contig_size:
            contigs_to_save.append({'contig_name': contig_name, 'contig_sequence': contig_sequence})
            current_size += contig_length
        else:
            i += 1
            contig_file_path = os.path.join('tmp', f'{name}_{i}_filtered.fasta')
            with open(contig_file_path, 'w') as contig_file:
                for saved_contig in contigs_to_save:
                    contig_file.write(f">{saved_contig['contig_name']}\n")
                    contig_file.write(f"{saved_contig['contig_sequence']}\n")
            contigs_to_save = []
            current_size = 0
    
    # Handle the last batch of contigs
    if contigs_to_save:
        i += 1
        contig_file_path = os.path.join('tmp', f'{name}_{i}_filtered.fasta')
        with open(contig_file_path, 'w') as contig_file:
            for saved_contig in contigs_to_save:
                contig_file.write(f">{saved_contig['contig_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")

def faa_into_chunk(faa_file, contig_size):
    contigs_to_save = []
    current_size = 0
    i = 0
    name = os.path.splitext(os.path.basename(str(cfg.genome_path)))[0]
    for record in pyfastx.Fasta(str(faa_file)):
        orf_name = str(record.name).split()[0]
        contig_sequence = str(record.seq)
        contig_length = len(record.seq)
        
        if current_size + contig_length <= contig_size:
            contigs_to_save.append({'orf_name': orf_name, 'contig_sequence': contig_sequence})
            current_size += contig_length
        else:
            i += 1
            contig_file_path = os.path.join('tmp', f'{name}_{i}_filtered.faa')
            with open(contig_file_path, 'w') as contig_file:
                for saved_contig in contigs_to_save:
                    contig_file.write(f">{saved_contig['orf_name']}\n")
                    contig_file.write(f"{saved_contig['contig_sequence']}\n")
            contigs_to_save = []
            current_size = 0
    
    # Handle the last batch of contigs
    if contigs_to_save:
        i += 1
        contig_file_path = os.path.join('tmp', f'{name}_{i}_filtered.faa')
        with open(contig_file_path, 'w') as contig_file:
            for saved_contig in contigs_to_save:
                contig_file.write(f">{saved_contig['orf_name']}\n")
                contig_file.write(f"{saved_contig['contig_sequence']}\n")


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
    check_plasmid_true(contig, 'is_circular')
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
                            'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                            contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                        )
                else:
                    hits_per_pos[hit_pos] = hit
                    log.info(
                        'inc-type: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    contig['inc_types'] = list(hits_per_pos.values())
    log.info('inc-type: contig=%s, # inc-types=%s', contig['id'], len(contig['inc_types']))
    check_plasmid_length(contig, 'inc_types')
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
                    'rRNAs: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
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
                    'ref plasmids: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
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
                    'oriT: hit! contig=%s, id=%s, c-start=%d, c-end=%d, coverage=%f, identity=%f',
                    contig['id'], hit['orit']['id'], hit['contig_start'], hit['contig_end'], hit['coverage'], hit['identity']
                )
    log.info('oriT: contig=%s, # oriT=%s', contig['id'], len(contig['orit_hits']))
    check_plasmid_length(contig, 'orit_hits')
    return


def filter_contig(contig):
    """Apply heuristic filters based on contig information."""
    if(contig['type'] == "plasmid"):
        return True
    elif(contig['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD):
        log.debug('filter: RDS > SPT! contig=%s', contig['id'])
        return True
    elif(contig['protein_score'] >= pc.RDS_CONSERVATIVE_THRESHOLD
            and len(contig['plasmid_hits']) > 0
            and len(contig['rrnas']) == 0):
        log.debug('filter: RDS > CT & plasmid hits & no rRNAs! contig=%s', contig['id'])
        return True
    else:
        return False
        
        
def filter_contig_meta(contig):
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
    if(contig['protein_score'] >= pc.RDS_SPECIFICITY_THRESHOLD):
        log.debug('filter: RDS > SPT! contig=%s', contig['id'])
        return True

    # include all contigs with mediocre protein scores but additional blast hit evidence without rRNAs
    if(contig['protein_score'] >= pc.RDS_CONSERVATIVE_THRESHOLD
            and len(contig['plasmid_hits']) > 0
            and len(contig['rrnas']) == 0):

        log.debug('filter: RDS > CT & plasmid hits & no rRNAs! contig=%s', contig['id'])
        return True

    return False
    

def train_gene_prediction(contigs): ##if(not pyrodigal_metamode):
    training_info = None
    log.info('create prodigal training info object: meta=%s, not pyrodigal_metamode')
    orf_finder = pyrodigal.OrfFinder(meta= False, closed=True)
    seqs = [contig['sequence'] for id, contig in contigs.items()]
    combined_sequence = ''.join(seqs)
    bytes_combined_sequence = combined_sequence.encode('utf-8')
    training_info = orf_finder.train(bytes_combined_sequence) ##  list comprehsion
    return training_info

def predict_orfs_py(contigs, record, proteins_path, training_info=None):
    """Predict open reading frames with Pyrodigal."""
    genome_size = sum([v['length'] for k, v in contigs.items()])
    pyrodigal_metamode = (cfg.module == 'metagenomic') or (cfg.characterize and genome_size < 50000) or (20000 > genome_size)
    log.info('ORFs: execute pyrodigal in meta mode! characterize=%s, genome-size=%d, metagenome=%s', cfg.characterize, genome_size, cfg.module)
    orf_finder = pyrodigal.OrfFinder(training_info, meta=pyrodigal_metamode, closed=True, mask=True)
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

def construct_hit_dict(hits, hit, orf, include_hmm_id=False):
    hit = {
        'type': hits.query_name.decode(),
        'start': int(orf['start']),
        'end': int(orf['end']),
        'strand': orf['strand'],
        'bitscore': hit.score,
        'evalue': hit.evalue
    }
    if include_hmm_id:
        hit['hmm-id'] = hit.hmm_name.decode()
    return hit
    

def search_amr_genes_py(contigs, filteredProteinsPath):
    """Search for AMR genes."""

    hmm_profile_path = str(cfg.db_path.joinpath('ncbifam-amr'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, bit_cutoffs="trusted"):
            for hit in hits:
                if hit.included:
                    contig_id = hit.name.decode()
                    tmp = contig_id.rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][int(tmp[1])]
                    hit = construct_hit_dict(hits, hit, orf, include_hmm_id= True)
                    contig['amr_hits'].append(hit)
                    log.info(
                        'AMRs: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('AMRs: contig=%s, # AMRs=%s', contig['id'], len(contig['amr_hits']))
    return


def search_replication_genes_py(contigs, filteredProteinsPath):
    """Search for replication genes."""

    hmm_profile_path = str(cfg.db_path.joinpath('replication'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-100):
            for hit in hits:
                if hit.included:
                    contig_id = hit.name.decode()
                    tmp = contig_id.rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][int(tmp[1])]
                    hit = construct_hit_dict(hits, hit, orf)
                    contig['replication_hits'].append(hit)
                    log.info(
                        'rep genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
                    check_plasmid_length(contig, 'replication_hits')
    log.info('rep genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['replication_hits']))
    return

def search_mobilization_genes_py(contigs, filteredProteinsPath):
    """Search for mobilization genes."""

    hmm_profile_path = str(cfg.db_path.joinpath('mobilization'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-10):
            for hit in hits:
                if hit.included:
                    contig_id = hit.name.decode()
                    tmp = contig_id.rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][int(tmp[1])]
                    hit = construct_hit_dict(hits, hit, orf)
                    contig['mobilization_hits'].append(hit)
                    log.info(
                        ' mob genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
                    check_plasmid_length(contig, 'mobilization_hits')
    log.info('mob genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['mobilization_hits']))
    return


def search_conjugation_genes_py(contigs, filteredProteinsPath):
    """Search for conjugation genes."""

    hmm_profile_path = str(cfg.db_path.joinpath('conjugation'))
    sequence_database_path = str(filteredProteinsPath)  # Replace with your sequence database file path
    with pyhmmer.easel.SequenceFile(sequence_database_path, digital=True) as seq_file:
        proteins = seq_file.read_block()
    with pyhmmer.plan7.HMMFile(hmm_profile_path) as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, proteins, cpus=1, E=1e-100):
            for hit in hits:
                if hit.included:
                    contig_id = hit.name.decode()
                    tmp = contig_id.rsplit('_', 1)
                    contig_id = tmp[0]
                    contig = contigs[contig_id]
                    orf = contig['orfs'][int(tmp[1])]
                    hit = construct_hit_dict(hits, hit, orf)
                    contig['conjugation_hits'].append(hit)
                    log.info(
                        'conj genes: hit! contig=%s, type=%s, start=%d, end=%d, strand=%s',
                        contig['id'], hit['type'], hit['start'], hit['end'], hit['strand']
                    )
    log.info('conj genes: contig=%s, # mob-genes=%s', contig['id'], len(contig['conjugation_hits']))
    return
