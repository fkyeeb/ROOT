import sys,os

args = sys.argv
if len(args) < 5:
    print("Usage: python sort_blast_result.py <blast_result_file> <genomic_seq2cluster.tsv_location> <all_or_best> <fasta_input_folder> <fasta_output_folder-optional> <tsv_output_folder-optional>")
    sys.exit(1)

def adjust_folder_name(folder):
    if folder[-1] != "/":
        folder += "/"
    return folder

blast_result_file = args[1]
seq2cluster_file = adjust_folder_name(args[2]) + "genomic_seq2cluster.tsv"
all_or_best = args[3]
fasta_input_folder = adjust_folder_name(args[4])

if len(args) > 5:
    fasta_output_folder = adjust_folder_name(args[5])
else:
    fasta_output_folder = ""

if len(args) > 6:
    tsv_output_folder = adjust_folder_name(args[6])
else:
    tsv_output_folder = ""

seq2cluster = {}
with open(seq2cluster_file,"r") as file:
    for line in file:
        splt = line.strip().split("\t")
        seq2cluster[splt[0]] = splt[1]

gene_cluster_seq = {}
with open(blast_result_file,"r") as file:
    for i in file:
        splt = i.strip().split("\t")
        genomic_seq = splt[1]
        qry_len = int(splt[3])
        sub_len = int(splt[4])
        aln_len = int(splt[7])
        sub_perc = aln_len / sub_len
        qry_perc = aln_len / qry_len
        gene = splt[0].split("-")[1]
        cluster = seq2cluster[genomic_seq]
        if gene not in gene_cluster_seq:
            gene_cluster_seq[gene] = {}
        if cluster not in gene_cluster_seq[gene]:
            gene_cluster_seq[gene][cluster] = {}
        if genomic_seq not in gene_cluster_seq[gene][cluster]:
            gene_cluster_seq[gene][cluster][genomic_seq] = [qry_perc, sub_perc]
        else:
            gene_cluster_seq[gene][cluster][genomic_seq][0] += qry_perc
            gene_cluster_seq[gene][cluster][genomic_seq][1] += sub_perc

# filter through results
for gene in list(gene_cluster_seq.keys())[:]:
    for cluster in list(gene_cluster_seq[gene].keys())[:]:
        rm = True
        for aln_perc in gene_cluster_seq[gene][cluster].values():
            if aln_perc[0] >= 0.9 or aln_perc[1] >= 0.5:
                rm = False
                break
        if rm:
            del gene_cluster_seq[gene][cluster]

# save results
with open(f"{tsv_output_folder}blast_result_filtered.tsv","w") as output:
    output.write("gene\tcluster\tgenomic_seq\tqry_perc\tsub_perc\n")
    for gene in gene_cluster_seq:
        for cluster in gene_cluster_seq[gene]:
            for seq, aln_perc in gene_cluster_seq[gene][cluster].items():
                output.write(gene + "\t" + cluster + "\t" + seq + "\t" + f"{aln_perc[0]:.2f}" + "\t" + f"{aln_perc[1]:.2f}" + "\n")

# group clusters into fasta files
for gene in gene_cluster_seq:
    with open(os.path.join(fasta_output_folder, gene + ".fa"),"w") as output:
        if all_or_best == "best":
            gene_cluster_seq[gene] = dict(sorted(gene_cluster_seq[gene].items(), key=lambda x: max(aln_perc[0] for aln_perc in x[1].values()), reverse=True)[:1])
        for cluster in gene_cluster_seq[gene]:
            with open(os.path.join(fasta_input_folder, cluster),"r") as f:
                out_str = f.read()
                output.write(out_str)