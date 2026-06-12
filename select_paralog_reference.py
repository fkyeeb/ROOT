import sys,os
import tree_reader as t

def adjust_folder_name(folder):
    if folder[-1] != "/":
        folder += "/"
    return folder

args = sys.argv
outgroup_file = args[1]
tree_folder = adjust_folder_name(args[2])
fasta_folder = adjust_folder_name(args[3])
gene_list_file = args[4]
fasta_file_ending = args[5]

if len(args) > 6:
    num_refs = int(args[6])
else:
    num_refs = 3

if len(args) > 7:
    ref_file_name = args[7]
else:
    ref_file_name = "hybpiper_paralog_reference.fa"

with open(outgroup_file,"r") as f:
    tip2ignore = set([line.strip() for line in f])

with open(gene_list_file,"r") as f:
    genes = [line.strip() for line in f]

tree_folder = tree_folder + "unpruned/"
potential_references = {}
for file in os.listdir(tree_folder):
    if file.endswith("rooted_processed_id.tre"):
        path = os.path.join(tree_folder, file)
        gene = file.split(".")[0]
        with open(path,"r") as f:
            tree = t.read_tree_string(f.readline().strip())
        potential_references[gene] = []
        for node in tree.iternodes():
            if node.label == "D":
                for child in node.children:
                    lvs = []
                    for tip in node.lvsnms():
                        taxon = tip.split("@")[0]
                        if taxon not in tip2ignore:
                            lvs.append(tip)
                    if len(lvs) > 0:
                        potential_references[gene].append(lvs)
        max_l = max([len(lvs) for lvs in potential_references[gene]], default=0)
        all_lvs = set(tree.lvsnms()) - tip2ignore
        if max_l < (len(all_lvs) / 2):
            potential_references[gene].append(all_lvs)
print(f"Number of trees read: {len(potential_references)}")
ref_str = ""
total_count = 0
total_genes = 0
for gene in genes:
    with open(os.path.join(fasta_folder, f"{gene}.{fasta_file_ending}"),"r") as f:
        seqs = {}
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq != "":
                    seqs[seq_name] = seq
                    seq = ""
                seq_name = line.strip()[1:]
            else:
                seq += line.strip()
        seqs[seq_name] = seq
    already_added = []
    count = 0
    # write longest seqs if no potential references
    if gene not in potential_references:
        seqs = dict(sorted(seqs.items(), key=lambda item: len(item[1]), reverse=True))
        for seq_name in list(seqs.keys())[:num_refs]:
            if not seq_name.endswith(f"-{gene}"):
                seq_name_long = f"{seq_name}-{gene}"
            out = f">{seq_name_long}\n{seqs[seq_name]}\n"
            ref_str += out
            count += 1
        total_genes += 1
    else:
        seq_lists = potential_references[gene]
        max_l = max([len(seq_list) for seq_list in seq_lists])
        for seq_list in seq_lists:
            seqs_subset = {taxon: seqs[taxon] for taxon in seq_list if taxon in seqs}
            seqs_subset = dict(sorted(seqs_subset.items(), key=lambda item: len(item[1]), reverse=True))
            if len(seqs_subset) < max_l:
                ref_num4subset = max(1, round(len(seqs_subset) / 10))
                temp_num_refs = min(num_refs, ref_num4subset)
            else:
                temp_num_refs = num_refs
            for seq_name in list(seqs_subset.keys())[:temp_num_refs]:
                if not seq_name.endswith(f"-{gene}"):
                    seq_name_long = f"{seq_name}-{gene}"
                if seq_name_long in already_added:
                    continue
                already_added.append(seq_name_long)
                out = f">{seq_name_long}\n{seqs_subset[seq_name]}\n"
                ref_str += out
                count += 1
        total_genes += 1
    total_count += count
    print(f"{gene}: {count} references")

# output refs
with open(ref_file_name,"w") as f:
    f.write(ref_str)
print(f"Total genes: {total_genes}")
print(f"Total references: {total_count}")
