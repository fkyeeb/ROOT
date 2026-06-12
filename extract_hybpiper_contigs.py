#!/usr/bin/env python

import sys
import re
import os
import math
from collections import Counter


def parse_allowed_contigs(filepath):
    seen = set()
    allowed = set()
    with open(filepath, 'r') as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 14:
                continue
            if 'query_id' in fields:
                continue
            query_id = fields[1].strip()
            hit_id   = fields[3].strip()
            if not query_id or not hit_id:
                continue
            if hit_id in ('Not performed', 'No Hits remaining', ''):
                continue
            dedup_key = (query_id, hit_id, fields[6].strip())
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            allowed.add(hit_id)
    return allowed


def kmer_complexity(seq, k=15):
    seq = seq.upper()
    n_kmers = len(seq) - k + 1
    if n_kmers <= 0:
        return 0.0
    counts = Counter(seq[i:i+k] for i in range(n_kmers))
    total = sum(counts.values())
    H = -sum((c / total) * math.log2(c / total) for c in counts.values())
    return H / math.log2(min(4**k, total))


def get_coverage(contig_id):
    m = re.search(r'_cov_([0-9]+(?:\.[0-9]+)?)', contig_id)
    return float(m.group(1)) if m else 0.0


def parse_exonerate_fasta(filepath, allowed_contigs, sample, gene):
    """
    Returns a list of dicts, one per passing contig:
        {
          'header':     '>sample-gene-contig_id',
          'sequence':   str,
          'contig_id':  str,
          'q_start_aa': int,
          'q_end_aa':   int,
          't_start':    int,
          't_end':      int,
        }
    Only the best-scoring alignment per contig is kept, then filtered by
    complexity and coverage.
    """
    best = {}  # contig_id -> dict
    current_score = None
    current_header = None
    current_seq_lines = []
    in_sequence = False

    score_re = re.compile(r'Raw score:\s*(\d+)')

    def _safe_int(s):
        s = s.strip()
        return int(s) if re.fullmatch(r'-?\d+', s) else 0

    def _commit(header, seq_lines, score):
        if header is None or score is None:
            return
        fields = header.split(',')
        if len(fields) != 8:
            return
        contig_id = fields[0]
        if contig_id not in allowed_contigs:
            return
        if contig_id not in best or score > best[contig_id]['score']:
            best[contig_id] = dict(
                score=score,
                sequence=''.join(seq_lines),
                contig_id=contig_id,
                q_start_aa=_safe_int(fields[2]),
                q_end_aa=_safe_int(fields[3]),
                t_start=_safe_int(fields[6]),
                t_end=_safe_int(fields[7]),
            )

    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            m = score_re.search(line)
            if m:
                if in_sequence and current_header is not None:
                    _commit(current_header, current_seq_lines, current_score)
                    in_sequence = False
                    current_header = None
                    current_seq_lines = []
                current_score = int(m.group(1))
                continue
            if line.startswith('>'):
                if in_sequence and current_header is not None:
                    _commit(current_header, current_seq_lines, current_score)
                current_header = line[1:]
                current_seq_lines = []
                in_sequence = True
                continue
            if in_sequence:
                stripped = line.strip()
                if stripped == '':
                    _commit(current_header, current_seq_lines, current_score)
                    in_sequence = False
                    current_header = None
                    current_seq_lines = []
                else:
                    current_seq_lines.append(stripped)

    if in_sequence and current_header is not None:
        _commit(current_header, current_seq_lines, current_score)

    if not best:
        return []

    # Filter: low-complexity
    best = {cid: rec for cid, rec in best.items()
            if kmer_complexity(rec['sequence'], k=15) > 0.99}
    if not best:
        return []

    # Filter: coverage < 10% of max
    max_cov = max(get_coverage(cid) for cid in best)
    min_cov = max_cov * 0.1
    best = {cid: rec for cid, rec in best.items()
            if get_coverage(cid) >= min_cov}

    result = []
    for cid, rec in best.items():
        result.append({
            'header':     f">{sample}-{gene}-{cid}",
            'sequence':   rec['sequence'],
            'contig_id':  cid,
            'q_start_aa': rec['q_start_aa'],
            'q_end_aa':   rec['q_end_aa'],
            't_start':    rec['t_start'],
            't_end':      rec['t_end'],
        })
    return result


def adjust_folder_name(folder):
    if folder[-1] != "/":
        folder += "/"
    return folder


def main():
    args = sys.argv
    sample_list_file       = args[1]
    gene_list_file         = args[2]
    genomic_seq_folder     = adjust_folder_name(args[3])
    hybpiper_output_folder = adjust_folder_name(args[4])
    output_folder          = adjust_folder_name(args[5])
    fasta_file_ending      = args[6] if len(args) > 6 else "fa"

    samples = [l.strip() for l in open(sample_list_file)]
    genes   = [l.strip() for l in open(gene_list_file)]

    positions_tsv_path = os.path.join(output_folder, "all_contig_positions.tsv")
    os.makedirs(output_folder, exist_ok=True)

    all_contigs_path = os.path.join(output_folder, "all_contigs.fa")

    with open(positions_tsv_path, 'w') as pos_out, \
         open(all_contigs_path, 'w') as all_out:

        pos_out.write("sample\tgene\tcontig_id\tq_start_aa\tq_end_aa\tt_start\tt_end\n")

        for sample in samples:
            sample_dir = os.path.join(output_folder, sample)
            os.makedirs(sample_dir, exist_ok=True)

            for gene in genes:
                exonerate_stats_file = (f"{hybpiper_output_folder}{sample}/{gene}"
                                        f"/{sample}/exonerate_stats.tsv")
                contig_file          = (f"{hybpiper_output_folder}{sample}/{gene}"
                                        f"/{sample}/exonerate_results.fasta")

                if not os.path.exists(exonerate_stats_file):
                    continue

                allowed_contigs = parse_allowed_contigs(exonerate_stats_file)
                contigs = parse_exonerate_fasta(contig_file, allowed_contigs, sample, gene)

                if not contigs:
                    continue

                print(f"Processing {sample} for {gene} ({len(contigs)} contig(s))")

                for rec in contigs:
                    all_out.write(f"{rec['header']}\n{rec['sequence']}\n")
                    pos_out.write(
                        f"{sample}\t{gene}\t{rec['contig_id']}\t"
                        f"{rec['q_start_aa']}\t{rec['q_end_aa']}\t"
                        f"{rec['t_start']}\t{rec['t_end']}\n"
                    )

                # Per-gene file (ref + contigs) for tree building, only if >1 contig
                if len(contigs) > 1:
                    gene_fa = os.path.join(sample_dir, f"{gene}.fa")
                    with open(gene_fa, 'w') as out:
                        with open(f"{genomic_seq_folder}{gene}.{fasta_file_ending}", 'r') as ref:
                            out_str = ref.read()
                            if not out_str.endswith('\n'):
                                out_str += '\n'
                            out.write(out_str)
                        for rec in contigs:
                            out.write(f"{rec['header']}\n{rec['sequence']}\n")

    print(f"Finished processing {len(samples)} samples")
    print(f"All contigs written to {all_contigs_path}")
    print(f"Positions written to {positions_tsv_path}")


if __name__ == '__main__':
    main()