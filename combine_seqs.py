import sys
import os

def adjust_folder(folder):
    if folder[-1] != "/":
        folder += "/"
    return folder

def main():
    args = sys.argv
    if len(args) < 6:
        print("Usage: python combine_seqs.py <sample_list_file> <gene_list_file> <ref_seq_folder> <contig_folder> <output_folder> [ref_fasta_file_ending]")
        sys.exit(1)
    
    sample_list_file = args[1]
    gene_list_file = args[2]
    ref_seq_dir  = adjust_folder(args[3])
    contig_dir   = adjust_folder(args[4])
    output_dir     = adjust_folder(args[5])
    ref_fasta_ending = args[6] if len(args) > 6 else "fa"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = [l.strip() for l in open(sample_list_file)]
    genes   = [l.strip() for l in open(gene_list_file)]

    for gene in genes:
        ref_fa = os.path.join(ref_seq_dir, f"{gene}.{ref_fasta_ending}")
        with open(output_dir + gene + ".fa", 'w') as out:
            with open(ref_fa, 'r') as ref:
                out_str = ref.read()
                if not out_str.endswith('\n'):
                    out_str += '\n'
                out.write(out_str)
            for sample in samples:
                contig_fa = os.path.join(contig_dir, sample, f"{gene}_stitched.fa")
                if os.path.exists(contig_fa):
                    with open(contig_fa, 'r') as contig:
                        out_str = contig.read()
                        if not out_str.endswith('\n'):
                            out_str += '\n'
                        out.write(out_str)

if __name__ == '__main__':
    main()