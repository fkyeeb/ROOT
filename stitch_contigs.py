#!/usr/bin/env python
"""
stitch_contigs.py

Stitch contig groups into consensus sequences using pre-computed group
assignments and reference mapping positions.

Usage:
    python stitch_contigs.py <sample_list> <gene_list> \
                             <all_contigs_dir> \
                             <contig_groups_tsv> \
                             <all_contig_positions_tsv> \
                             <output_folder>

Overlap handling
----------------
Contigs within a group are sorted by q_start_aa.  If contig B's q_start_aa
falls before contig A's q_end_aa, they overlap on the reference.

  Case 1 – B is fully contained inside A (B.q_end <= A.q_end):
      B is dropped entirely; it adds no new reference coverage.

  Case 2 – B partially overlaps A's right end (B.q_start < A.q_end <= B.q_end):
      The overlap in nucleotides is (A.q_end - B.q_start) * 3.
      The higher-coverage contig keeps the overlapping bases; the lower-coverage
      one is trimmed at that boundary:
        - if A has higher coverage: B's sequence is left-trimmed by overlap_nt bases
        - if B has higher coverage: A's sequence (already appended) is right-trimmed
          by overlap_nt bases before appending B

Both cases are recorded in stitched_info.tsv with an "overlap_action" column.

Output per sample
-----------------
  <output_folder>/<sample>/<gene>_stitched.fa
  <output_folder>/<sample>/stitched_all.fa
  <output_folder>/<sample>/stitched_info.tsv
"""

import sys
import os
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_fasta(filepath):
    """Return {header_without_gt: sequence}."""
    seqs = {}
    header = None
    buf = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(buf)
                header = line[1:]
                buf = []
            else:
                buf.append(line.strip())
    if header is not None:
        seqs[header] = ''.join(buf)
    return seqs


def load_positions(positions_tsv):
    """Return {(sample, gene, contig_id): (q_start_aa, q_end_aa, t_start, t_end)}."""
    pos = {}
    with open(positions_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('sample') or not line:
                continue
            fields = line.split('\t')
            if len(fields) < 7:
                continue
            sample, gene, contig_id = fields[0], fields[1], fields[2]
            try:
                q_start, q_end = int(fields[3]), int(fields[4])
                t_start, t_end = int(fields[5]), int(fields[6])
            except ValueError:
                q_start = q_end = t_start = t_end = 0
            pos[(sample, gene, contig_id)] = (q_start, q_end, t_start, t_end)
    return pos


def load_groups(groups_tsv):
    """Return {(sample, gene): [[contig_id, ...], ...]}."""
    groups = defaultdict(list)
    with open(groups_tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            sample, gene = fields[0], fields[1]
            groups[(sample, gene)].append(fields[2].split(','))
    return groups


def get_coverage(contig_id):
    """Extract cov_XXX from a SPAdes contig name."""
    m = re.search(r'_cov_([0-9]+(?:\.[0-9]+)?)', contig_id)
    return float(m.group(1)) if m else 0.0


# ---------------------------------------------------------------------------
# Stitching — returns (stitched_sequence, list_of_info_dicts)
# ---------------------------------------------------------------------------

def stitch_group(contig_seqs, group, positions, sample, gene):
    """
    Stitch contigs in `group` ordered by q_start_aa.
    Handles gaps (N-padding) and overlaps (trim lower-coverage side).

    Returns:
        stitched_seq  : str
        info_rows     : list of dicts, one per contig considered, with keys:
                          contig_id, action, overlap_action,
                          stitched_base_start, stitched_base_end,
                          q_start_aa, q_end_aa, coverage
    """
    def get_seq(cid):
        return contig_seqs.get(f"{sample}-{gene}-{cid}", "")

    def get_pos(cid):
        return positions.get((sample, gene, cid), (0, 0, 0, 0))

    ordered = sorted(group, key=lambda c: get_pos(c)[0])

    pieces       = []    # sequence strings to join at the end
    cursor       = 0     # current total length of joined pieces
    info_rows    = []    # one dict per contig (including dropped)
    prev_end_aa  = None
    prev_cid     = None
    prev_row_idx = None  # index into info_rows for the last used contig

    for cid in ordered:
        seq            = get_seq(cid)
        q_start, q_end = get_pos(cid)[:2]
        cov            = get_coverage(cid)

        if not seq:
            info_rows.append(dict(
                contig_id=cid, action="dropped",
                overlap_action="no_sequence",
                stitched_base_start="", stitched_base_end="",
                q_start_aa=q_start, q_end_aa=q_end, coverage=cov,
            ))
            continue

        overlap_action_current = "none"

        if prev_end_aa is not None and q_start < prev_end_aa:
            overlap_nt = (prev_end_aa - q_start) * 3

            if q_end <= prev_end_aa:
                # Case 1: fully contained — drop current
                info_rows.append(dict(
                    contig_id=cid, action="dropped",
                    overlap_action="fully_contained_in_previous",
                    stitched_base_start="", stitched_base_end="",
                    q_start_aa=q_start, q_end_aa=q_end, coverage=cov,
                ))
                continue

            # Case 2: partial overlap — higher coverage wins the disputed bases
            prev_cov = get_coverage(prev_cid)
            if prev_cov >= cov:
                # Previous wins: left-trim the current contig
                seq = seq[overlap_nt:]
                overlap_action_current = f"left_trimmed_{overlap_nt}nt_lower_coverage"
            else:
                # Current wins: right-trim the previous contig's piece
                pieces[-1] = pieces[-1][:-overlap_nt] if overlap_nt < len(pieces[-1]) else ""
                cursor -= overlap_nt
                # Update the previous contig's recorded end in its info row
                prev_row = info_rows[prev_row_idx]
                prev_row["stitched_base_end"]  = prev_row["stitched_base_end"] - overlap_nt
                prev_row["overlap_action"]     = f"right_trimmed_{overlap_nt}nt_lower_coverage"
                overlap_action_current = "none"

        # (gaps between contigs are not padded — contigs are concatenated directly)

        # Append current contig
        base_start = cursor + 1        # 1-based
        base_end   = cursor + len(seq)
        pieces.append(seq)
        cursor += len(seq)

        prev_row_idx = len(info_rows)
        info_rows.append(dict(
            contig_id=cid, action="used",
            overlap_action=overlap_action_current,
            stitched_base_start=base_start,
            stitched_base_end=base_end,
            q_start_aa=q_start, q_end_aa=q_end, coverage=cov,
        ))

        prev_end_aa = max(prev_end_aa or 0, q_end)
        prev_cid    = cid

    return "".join(pieces), info_rows


def get_longest(contig_seqs, contig_ids, sample, gene):
    return max(contig_ids,
               key=lambda c: len(contig_seqs.get(f"{sample}-{gene}-{c}", '')))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

INFO_HEADER = (
    "output_seq\tgene\tcontig_id\taction\toverlap_action\t"
    "stitched_base_start\tstitched_base_end\t"
    "q_start_aa\tq_end_aa\tcoverage\n"
)


def adjust_folder(folder):
    return folder if folder.endswith('/') else folder + '/'


def main():
    args = sys.argv
    if len(args) < 6:
        print(__doc__)
        sys.exit(1)

    sample_list_file = args[1]
    gene_list_file   = args[2]
    all_contigs_dir  = adjust_folder(args[3])
    groups_tsv       = args[4]
    output_folder    = adjust_folder(args[5])
    max_copies       = int(args[6]) if len(args) > 6 else None   # None = no limit
    min_overlap_bp   = int(args[7]) if len(args) > 7 else 50

    samples = [l.strip() for l in open(sample_list_file)]
    genes   = [l.strip() for l in open(gene_list_file)]

    groups    = load_groups(groups_tsv)
    positions = load_positions(os.path.join(all_contigs_dir, "all_contig_positions.tsv"))

    all_contigs_fa = os.path.join(all_contigs_dir, "all_contigs.fa")
    if not os.path.exists(all_contigs_fa):
        print(f"[ERROR] all_contigs.fa not found at {all_contigs_fa}")
        sys.exit(1)

    print("Loading all_contigs.fa ...")
    all_seqs = read_fasta(all_contigs_fa)

    # Index contig_ids per (sample, gene)
    gene_contig_ids = defaultdict(list)
    for header in all_seqs:
        parts = header.split('-', 2)
        if len(parts) == 3:
            gene_contig_ids[(parts[0], parts[1])].append(parts[2])

    # copy_counts[gene][sample] = number of stitched copies output
    copy_counts = {gene: {sample: 0 for sample in samples} for gene in genes}

    for sample in samples:
        out_sample_dir = os.path.join(output_folder, sample)
        os.makedirs(out_sample_dir, exist_ok=True)

        stitched_all_records = []
        all_info_rows = []          # (output_seq_name, gene, info_dict)

        for gene in genes:
            contig_ids  = gene_contig_ids.get((sample, gene))
            if not contig_ids:
                continue

            gene_groups = groups.get((sample, gene))
            gene_records = []       # [(out_header, seq), ...]

            if gene_groups:
                # Stitch all groups first
                stitched = []   # [(group, seq, info_rows), ...]
                for group in gene_groups:
                    seq, info_rows = stitch_group(all_seqs, group, positions, sample, gene)
                    if seq:
                        stitched.append((group, seq, info_rows))

                # Filter groups by reference overlap when there are multiple
                if len(stitched) > 1:
                    # Compute covered reference intervals for each group,
                    # using only contigs that were actually used (not dropped).
                    # Each interval is in nt (q coords * 3).
                    def group_covered_intervals(info_rows):
                        intervals = []
                        for row in info_rows:
                            if row["action"] == "used":
                                intervals.append(
                                    (row["q_start_aa"] * 3, row["q_end_aa"] * 3)
                                )
                        return intervals

                    covered = [group_covered_intervals(info_rows)
                               for _, _, info_rows in stitched]

                    def ref_overlap_bp(intervals_a, intervals_b):
                        """Sum of actual intersections across all interval pairs."""
                        total = 0
                        for a0, a1 in intervals_a:
                            for b0, b1 in intervals_b:
                                total += max(0, min(a1, b1) - max(a0, b0))
                        return total

                    n = len(stitched)
                    keep = set()
                    for i in range(n):
                        for j in range(i + 1, n):
                            if ref_overlap_bp(covered[i], covered[j]) >= min_overlap_bp:
                                keep.add(i)
                                keep.add(j)

                    if keep:
                        stitched = [stitched[i] for i in sorted(keep)]
                    else:
                        # No pair overlaps enough — keep only the largest group
                        largest = max(range(n), key=lambda i: len(stitched[i][0]))
                        stitched = [stitched[largest]]

                # Rank by sequence length (excluding Ns) and cap at max_copies
                stitched.sort(key=lambda x: len(x[1].replace("N", "")), reverse=True)
                if max_copies is not None:
                    stitched = stitched[:max_copies]

                n_kept = len(stitched)
                for i, (group, seq, info_rows) in enumerate(stitched, 1):
                    out_name   = f"{sample}@{i}-{gene}" if n_kept > 1 else f"{sample}@{gene}"
                    out_header = f">{out_name}"
                    gene_records.append((out_header, seq))
                    for row in info_rows:
                        all_info_rows.append((out_name, gene, row))

            else:
                # No group — longest single contig
                longest  = get_longest(all_seqs, contig_ids, sample, gene)
                seq      = all_seqs[f"{sample}-{gene}-{longest}"]
                out_name = f"{sample}@{gene}"
                gene_records.append((f">{out_name}", seq))
                cov = get_coverage(longest)
                q_start, q_end = positions.get((sample, gene, longest), (0, 0, 0, 0))[:2]
                all_info_rows.append((out_name, gene, dict(
                    contig_id=longest, action='used_as_singleton',
                    overlap_action='none',
                    stitched_base_start=1, stitched_base_end=len(seq),
                    q_start_aa=q_start, q_end_aa=q_end, coverage=cov,
                )))

            # Per-gene stitched FASTA
            gene_out = os.path.join(out_sample_dir, f"{gene}_stitched.fa")
            with open(gene_out, 'w') as fh:
                for out_header, seq in gene_records:
                    fh.write(f"{out_header}\n{seq}\n")

            copy_counts[gene][sample] = len(gene_records)
            stitched_all_records.extend(gene_records)
            copies = ', '.join(h.lstrip('>') for h, _ in gene_records)
            print(f"  {gene}: {copies}")

        # stitched_all.fa
        combined_out = os.path.join(out_sample_dir, "stitched_all.fa")
        with open(combined_out, 'w') as fh:
            for out_header, seq in stitched_all_records:
                fh.write(f"{out_header}\n{seq}\n")

        # stitched_info.tsv
        info_out = os.path.join(out_sample_dir, "stitched_info.tsv")
        with open(info_out, 'w') as fh:
            fh.write(INFO_HEADER)
            for out_name, gene, row in all_info_rows:
                fh.write(
                    f"{out_name}\t{gene}\t{row['contig_id']}\t"
                    f"{row['action']}\t{row['overlap_action']}\t"
                    f"{row['stitched_base_start']}\t{row['stitched_base_end']}\t"
                    f"{row['q_start_aa']}\t{row['q_end_aa']}\t{row['coverage']}\n"
                )

        print(f"[DONE] {sample}: {len(stitched_all_records)} records → stitched_all.fa + stitched_info.tsv")

    # Write copy-count summary TSV
    summary_path = os.path.join(output_folder, "copy_counts.tsv")
    with open(summary_path, "w") as fh:
        sep = "\t"
        nl  = "\n"
        fh.write(sep + sep.join(samples) + nl)
        totals = [sum(1 for gene in genes if copy_counts[gene][s] > 1) for s in samples]
        fh.write("total_genes_>1_copy" + sep + sep.join(map(str, totals)) + nl)
        for gene in genes:
            counts = [copy_counts[gene][s] for s in samples]
            fh.write(gene + sep + sep.join(map(str, counts)) + nl)
    print(f"Copy counts written to {summary_path}")
    print(f"Finished {len(samples)} samples.")


if __name__ == '__main__':
    main()