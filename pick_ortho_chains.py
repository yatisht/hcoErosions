
import argparse
import numpy as np
import os.path
import os
import math
import subprocess
import pickle
import sys
import tempfile
import gzip

CHAIN_DIR = 'chains/' 
CHAIN_ID_DIR = 'chain_ids/'
REF_GENE_DIR = 'ref_genes/'
BED12_TO_BED6 = 'bin/bedToExons'

def get_gene_chr_coords (bed12filename):
    genes_chr = {}
    genes_coords = {}

    with open(bed12filename, "r") as f:
        for line in f:
            words = line.split()
            gene_name = words[3] 
            gene_chr = words[0]
            genes_chr[gene_name] = gene_chr
    
    tmpBed = tempfile.NamedTemporaryFile(suffix='.bed')
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call([BED12_TO_BED6, "-cdsOnly", 
                               "-i", bed12filename, "-o", tmpBed.name]) 

    for line in tmpBed:
        words = line.split()
        gene_name = words[3] 
        gene_chr = words[0]
        start = int(words[1])
        end = int(words[2])
        coords = genes_coords.get(gene_name, [])
        coords.append([start, end])
        genes_coords[gene_name] = coords

    tmpBed.close()
    return (genes_chr, genes_coords)

def find_query_coords(chain_id, t_start, t_end, strand):
    chain_start = chains[chain_id][1]
    chain_end = chains[chain_id][2]
    if (chain_end < t_start) or (chain_start > t_end):
        return ['', 0, 0, 0]
    else:
        q_chr = chains[chain_id][4]
        t_curr = chain_start
        q_chr_len = chains[chain_id][5]
        q_curr = chains[chain_id][6]
        q_start = 0
        q_end = 0
        overlap_start = 0
        overlap_end = 0
        overlap = 0
        gaps_len = len(chain_t_gaps[chain_id])
        for (k, block) in enumerate(chain_blocks[chain_id]):
            overlap += max(0, min(t_curr + block, t_end) - max(t_curr, t_start))
            if (q_start == 0) and (t_curr + block >= t_start):
                q_start = q_curr + max(0, (t_start - t_curr))
            if (q_end == 0) and (t_curr + block >= t_end):
                q_end = q_curr + max(0, t_end - t_curr)
            t_curr += block
            q_curr += block
            if (k < gaps_len):
                t_curr += chain_t_gaps[chain_id][k]
                q_curr += chain_q_gaps[chain_id][k]
        if (q_end == 0):
            q_end = q_curr
        if strand == '+':
            return [q_chr, q_start, q_end, overlap]
        else:
            return [q_chr, q_chr_len - q_end, q_chr_len - q_start, overlap]

def find_orthologous_gene(t_chr, t_exon_starts, t_exon_ends, gene_id):
    gene_len = sum([(t_exon_ends[k] - t_exon_starts[k]) for k in \
                    range(len(t_exon_starts))])
    best_chain_q_chr = ''
    best_chain_q_strand = ''
    best_chain_q_starts = [0] * len(t_exon_starts)
    best_chain_q_ends = [0] * len(t_exon_ends)
    exon_avail = [1] * len(t_exon_ends)
    best_chain_len = 1
    level = 1
    best_chain_level = 1
    highest_overlap = 1
    second_highest_overlap = 0
    highest_chain_score = 1
    second_highest_chain_score = 0
    overlapping_chains = []
    best_chain_id = -1
    second_best_chain_id = -1
    highest_overlap_score = 0
    level1_score = 1
    level1_chr = ''
    level1_overlap = 0
    level1_chain_id = -1
    level1_gene_in_synteny = 1
    level1_q_starts = [0] * len(t_exon_starts)
    level1_q_ends = [0] * len(t_exon_ends)
    overlapping_exon_highest_chain_score = [0] * len(t_exon_starts)
    coding_start = t_exon_starts[0]
    coding_end = t_exon_ends[-1]
    if (coding_end <= coding_start):
        coding_start = t_exon_starts[-1]
        coding_end = t_exon_ends[0]
    chain_ids = []
    if t_chr in chain_intervals.keys():
        interval_start = (t_exon_starts[0] / interval_len)
        interval_end = 1 + (t_exon_ends[-1] / interval_len)
        if (interval_end <= interval_start):
            interval_start = (t_exon_starts[-1] / interval_len)
            interval_end = 1 + (t_exon_ends[0] / interval_len)
        for i in range(interval_start, interval_end):
            if (i < len(chain_intervals[t_chr])):
                chain_ids.extend(chain_intervals[t_chr][i])
    chains_keys = sorted(set(chain_ids))
    for chain_id in chains_keys:
        chain = chains[chain_id]
        chain_score = chain[-1]
        if (chain[0] == t_chr):
            chain_start = chain[1]
            chain_end = chain[2]
            chain_strand = chain[3]
            chain_overlap = min(coding_end, chain_end) - max(chain_start, \
                                                             coding_start)
            if (chain_overlap > 0):
                q_chr = chain[5]
                q_starts = []
                q_ends = []
                overlap_score = []
                overlaps = []
                for (k, exon_start) in enumerate(t_exon_starts):
                    exon_end = t_exon_ends[k]
                    coords = find_query_coords(chain_id, exon_start, exon_end, \
                                               chain_strand)
                    overlap = coords[3]
                    if (exon_avail[k] > 0):
                        if (coords[3] > (t_exon_ends[k] - t_exon_starts[k])/4):
                            exon_avail[k] = 0
                    else:
                        overlap = 0
                    q_starts.append(coords[1])
                    q_ends.append(coords[2])
                    overlap_score.append(overlap)
                    overlaps.append(coords[3])
                total_overlap_score = sum(overlap_score)
                total_overlap = sum(overlaps)
                overlapping_chains.append([chain_id, total_overlap, \
                                           chain_score])
                if (total_overlap_score > highest_overlap_score):
                    highest_overlap_score = total_overlap_score
                    best_chain_id = chain_id
                    best_chain_q_chr = chain[4]
                    best_chain_q_strand = chain[3]
                    best_chain_q_starts = q_starts
                    best_chain_q_ends = q_ends
                    best_chain_level = level
                    highest_overlap = total_overlap
                    highest_chain_score = chain_score
                    best_chain_len = chain_end - chain_start
                if (chain_overlap == coding_end - coding_start):
                    if (level == 1):
                        level1_score = chain_score
                        level1_overlap = float(total_overlap) / gene_len
                        level1_chain_id = chain_id
                        level1_chr = chain[4]
                        level1_q_starts = q_starts
                        level1_q_ends = q_ends
                    level = level+1
    for overlapping_chain in overlapping_chains:
        chain_id = overlapping_chain[0]
        total_overlap = overlapping_chain[1]
        chain_score = overlapping_chain[2]
        if ((2*total_overlap > highest_overlap) and \
            (chain_id != best_chain_id) and \
            (chain_score > second_highest_chain_score)):
            second_highest_chain_score = chain_score
            second_highest_overlap = total_overlap
            second_best_chain_id = chain_id
    gene_in_synteny = float(gene_len) / chain_lens.get(best_chain_id, 1)
    level1_gene_in_synteny = float(gene_len) / \
            chain_lens.get(level1_chain_id, 1)
    second_best_score_ratio = (float(second_highest_chain_score) / \
                               highest_chain_score) 
    second_best_overlap_ratio = (float(second_highest_overlap) / \
                                 highest_overlap) 
    level1_score_ratio = float(highest_chain_score) / level1_score
    return (best_chain_id, second_best_chain_id, \
            best_chain_q_chr, best_chain_q_starts, best_chain_q_ends, \
            best_chain_level, gene_in_synteny, second_best_score_ratio, \
            second_best_overlap_ratio, best_chain_q_strand, level1_chain_id, \
            level1_chr, level1_q_starts, level1_q_ends, \
            level1_gene_in_synteny, level1_score_ratio, level1_overlap)

def is_assigned(best_chain_id, best_chain_q_chr, best_chain_q_starts, \
                best_chain_q_ends, assigned):
    overlap = 0
    for coords in assigned[best_chain_q_chr]:
        for (k, q_start) in enumerate(best_chain_q_starts):
            q_end = best_chain_q_ends[k]
            overlap = max(0, min(coords[1], q_end) - max(coords[0], q_start))
            if ((overlap > 0) and (best_chain_id != coords[2])):
                return True
    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Picks orthologous chain from '
            'reference-query alignments for each gene in the reference gene \
                                     set.')
    parser.add_argument("-reference", type=str,
            help="reference assembly name")
    parser.add_argument("-query", type=str,
            help="query assembly name")
    parser.add_argument("-genes", type=str,
            help="reference gene set (bed12 file)")
    parser.add_argument("-geneTranscriptIds", type=str,
            help="reference gene set complete trancript ids")    
    parser.add_argument("-geneCanonicalTranscriptIds", type=str,
            help="reference gene set canonical trancript ids")    
    
    if len(sys.argv) <= 10:
        parser.print_help()
        sys.exit(1)
    
    args = vars(parser.parse_args())
    
    reference = args['reference']
    query = args['query']
    bed12filename = args['genes']
    gene_transcript_filename = args['geneTranscriptIds']
    gene_canonical_transcript_filename = args['geneCanonicalTranscriptIds']

    gene_transcript_dict = {}
    with open(gene_transcript_filename, 'r') as f:
        for line in f:
            words = line.split()
            l = gene_transcript_dict.get(words[0], [])
            l.append(words[1])
            gene_transcript_dict[words[0]] = l
    
    gene_canonical_transcript_dict = {}
    with open(gene_canonical_transcript_filename, 'r') as f:
        for line in f:
            words = line.split()
            gene_canonical_transcript_dict[words[0]] = words[1]
    
    synteny_thresh = 0.05
    second_best_thresh = 0.05

    overlap_thresh = 0.3
    score_thresh = 0.1
    
    min_chain_score = 30000
    interval_len = 10000
    
    chainfile = CHAIN_DIR + reference + '.' + query + '.all.chain.gz'
    
    chains = {} 
    chain_blocks = {} 
    chain_t_gaps = {} 
    chain_q_gaps = {} 
    chain_id = 0
    chain_intervals = {}
    
    assigned = {}
    not_assigned = []
    selected_chains = set([])
    
    with gzip.open(chainfile, "r") as f:
        for line in f:
            words = line.split()
            if (line[0:5] == 'chain'):
                chain_chr = words[2]
                if chain_chr not in chain_intervals.keys():
                    chain_intervals[chain_chr] = [[]] * (1 + (int(words[3]) / \
                                                              interval_len))
                chain_start = int(words[5])
                chain_end = int(words[6])
                chain_q_chr = words[7]
                chain_q_chr_len = int(words[8])
                if chain_q_chr not in assigned.keys():
                    assigned[chain_q_chr] = []
                chain_q_strand = words[9]
                chain_q_start = int(words[10])
                chain_q_end = int(words[11])
                chain_score = int(words[1])
                chain_id = int(words[12])
                if chain_score < min_chain_score:
                    break
                chains[chain_id] = [chain_chr, chain_start, chain_end, \
                                chain_q_strand, chain_q_chr, chain_q_chr_len, \
                                chain_q_start, chain_q_end, chain_score]
                chain_blocks[chain_id] = []
                chain_t_gaps[chain_id] = []
                chain_q_gaps[chain_id] = []
                interval_start = (chain_start / interval_len)
                interval_end = 1 + (chain_end / interval_len)
                for i in range(interval_start,   interval_end):
                    chain_intervals[chain_chr][i].append(chain_id)
            elif (chain_id > 0):
                if (len(words) > 0):
                    block_length = int(words[0])
                    chain_blocks[chain_id].append(block_length)
                if (len(words) > 1):
                    t_gap_length = int(words[1])
                    q_gap_length = int(words[2])
                    chain_t_gaps[chain_id].append(t_gap_length)
                    chain_q_gaps[chain_id].append(q_gap_length)

    chain_lens = {}
    for chain_id in chain_blocks.keys():
        chain_lens[chain_id] = sum(chain_blocks[chain_id])

    (genes_chr, genes_coords) = get_gene_chr_coords (bed12filename)
    transcript_chain_id = {}

    count = 0

    for (num, gene_id) in enumerate(gene_canonical_transcript_dict.keys()):
         if (num % 1000 == 0):
             print '- Processed ', num, ' genes.'
         transcript_id = gene_canonical_transcript_dict[gene_id]
         chr =  genes_chr[transcript_id]
         exon_starts = [genes_coords[transcript_id][k][0] for k in \
                        range(len(genes_coords[transcript_id]))]
         exon_ends = [genes_coords[transcript_id][k][1] for k in \
                      range(len(genes_coords[transcript_id]))]
         (best_chain_id, second_best_chain_id, best_chain_q_chr, \
          best_chain_q_starts, best_chain_q_ends, best_chain_level, \
          gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, \
          best_chain_q_strand, level1_chain_id, level1_chr, level1_q_starts, \
          level1_q_ends, level1_gene_in_synteny, level1_score_ratio, \
          level1_overlap) = (find_orthologous_gene(chr, exon_starts, \
                                                   exon_ends, transcript_id))
         if ((gene_in_synteny < synteny_thresh) and (second_best_score_ratio < \
                           second_best_thresh) and (best_chain_level == 1)):
             for tid in gene_transcript_dict[gene_id]:
                 transcript_chain_id[tid] = [str(best_chain_id)]
             selected_chains.update(set([best_chain_id]))
             for (k, q_start) in enumerate(best_chain_q_starts):
                 q_end = best_chain_q_ends[k]
                 if (q_end > q_start):
                     assigned[best_chain_q_chr].append((q_start, q_end, \
                                                        best_chain_id))
                 else: 
                     assigned[best_chain_q_chr].append((q_end, q_start, \
                                                        best_chain_id))
         else:
            not_assigned.append(gene_id)
    
    for gene_id in not_assigned:
         transcript_id = gene_canonical_transcript_dict[gene_id]
         chr =  genes_chr[transcript_id]
         exon_starts = [genes_coords[transcript_id][k][0] for k in \
                        range(len(genes_coords[transcript_id]))]
         exon_ends = [genes_coords[transcript_id][k][1] for k in \
                      range(len(genes_coords[transcript_id]))]
         (best_chain_id, second_best_chain_id, best_chain_q_chr, \
          best_chain_q_starts, best_chain_q_ends, best_chain_level, \
          gene_in_synteny, second_best_score_ratio, second_best_overlap_ratio, \
          best_chain_q_strand, level1_chain_id, level1_chr, level1_q_starts, \
          level1_q_ends, level1_gene_in_synteny, level1_score_ratio, \
          level1_overlap) = (find_orthologous_gene(chr, exon_starts, \
                                                   exon_ends, transcript_id)) 
         if (((gene_in_synteny < synteny_thresh) or ((best_chain_level==1) \
                    and (best_chain_id in selected_chains))) and \
             (second_best_score_ratio < second_best_thresh) and \
             (level1_score_ratio > score_thresh or level1_overlap <= \
              overlap_thresh)):
             if not(is_assigned(best_chain_id, best_chain_q_chr, \
                            best_chain_q_starts, best_chain_q_ends, assigned)):
                    selected_chains.update(set([best_chain_id]))
                    for tid in gene_transcript_dict[gene_id]:
                        transcript_chain_id[tid] = [str(best_chain_id)]
                    for (k, q_start) in enumerate(best_chain_q_starts):
                        q_end = best_chain_q_ends[k]
                        if (q_end > q_start):
                            assigned[best_chain_q_chr].append((q_start, q_end, \
                                                               best_chain_id))
                        else: 
                            assigned[best_chain_q_chr].append((q_end, q_start, \
                                                               best_chain_id))
         elif ((level1_overlap > overlap_thresh) and (level1_gene_in_synteny \
                                                      < synteny_thresh)):
             if not(is_assigned(level1_chain_id, level1_chr, level1_q_starts, \
                                level1_q_ends, assigned)):
                    selected_chains.update(set([level1_chain_id]))
                    for tid in gene_transcript_dict[gene_id]:
                        transcript_chain_id[tid] = [str(level1_chain_id)]
                    for (k, q_start) in enumerate(level1_q_starts):
                        q_end = level1_q_ends[k]
                        if (q_end > q_start):
                            assigned[level1_chr].append((q_start, q_end, \
                                                         level1_chain_id))
                        else: 
                            assigned[level1_chr].append((q_end, q_start, \
                                                         level1_chain_id))

    outfile = CHAIN_ID_DIR + query + '.p'
    pickle.dump(transcript_chain_id, open(outfile, 'w'))

