import os,sys
from collections import Counter
import tree_reader as t
from itertools import combinations

def adjust_folder_name(folder):
    if folder[-1] != "/":
        folder += "/"
    return folder

def check_sister(node):
    # check if the sister clade has any tips that are not contigs (i.e. reference sequences)
    if node.parent is None:
        return False
    for child in node.parent.children:
        if child != node:
            for tip in child.iternodes():
                if tip.istip and "length" not in tip.label:
                    return True
    return False

def check_if_all_contigs(node):
    for tip in node.iternodes():
        if tip.istip and "length" not in tip.label:
            return False
    return True

def trace(node, other_trace = None, other_l = None, contig_name = None):
    trace = {}
    depth = 0
    l = []
    contig_len = node.length
    while node is not None:
        if other_trace:
            if node in other_trace:
                num_nodes = depth + other_trace[node] - 1
                if num_nodes < 0:
                    num_nodes = 0
                # not double counting mrca
                total_len = sum(l[1:depth]) + sum(other_l[1:other_trace[node]])
                # not including mrca's branch length and the tips themselves
                if check_if_all_contigs(node):
                    if node.length >= contig_len * 8 or node.length >= 0.1:
                        return None, None
                return num_nodes, total_len
        if check_sister(node):
            trace[node] = depth
            l.append(node.length)
            depth += 1
        else:
            l.append(0)
        node = node.parent

    if other_trace:
        num_nodes = depth + len(other_trace)
        if num_nodes < 0:
            num_nodes = 0
        # not double counting mrca
        total_len = sum(l[1:depth]) + sum(other_l[1:])
        # not including mrca's branch length and the tips themselves
        return num_nodes, total_len
    
    return trace, l

def find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]  # path compression
        x = parent[x]
    return x

def union(parent, x, y):
    parent[find(parent, x)] = find(parent, y)

def main():
    args = sys.argv
    folder = adjust_folder_name(args[1])
    tree_file_ending = args[2]
    sample_list_file = args[3]
    if len(args) > 4:
        num_nodes_threshold = int(args[4])
    else:
        num_nodes_threshold = 5
    if len(args) > 5:
        total_len_threshold = float(args[5])
    else:
        total_len_threshold = 0.05

    samples = [sample for sample in open(sample_list_file).read().splitlines()]

    # out = open("contig_groups.tsv", "w")
    out = open(f"contig_groups-{args[1]}.tsv", "w")

    for sample in samples:
        sample_folder = folder + sample + "/"
        for file in os.listdir(sample_folder):
            if file.endswith(tree_file_ending):
                path = os.path.join(sample_folder, file)
                gene = file.strip(tree_file_ending)
                with open(path, 'r') as f:
                    tree = t.read_tree_string(f.readline().strip())
                contigs = []
                for tip in tree.iternodes():
                    if tip.istip and "length" in tip.label:
                        contigs.append(tip)
                if len(contigs) <= 1:
                    continue
                contig_names = [c.label.split("-")[-1] for c in contigs]
                parent = {name: name for name in contig_names}
                for contig1, contig2 in combinations(contigs, 2):
                    trace1, l1 = trace(contig1)
                    num_nodes, total_len = trace(contig2, trace1, l1, contig1.label.split("-")[-1])
                    contig1_name = contig1.label.split("-")[-1]
                    contig2_name = contig2.label.split("-")[-1]
                    # if gene == "6848" and sample == "HEPU_short":
                    #     print(f"{gene} {sample} {contig1.label} vs {contig2.label}: num_nodes={num_nodes}, total_len={total_len}")
                    if num_nodes is not None and total_len is not None:
                        if num_nodes <= num_nodes_threshold or total_len <= total_len_threshold:
                            union(parent, contig1_name, contig2_name)
                groups = {}
                for name in contig_names:
                    root = find(parent, name)
                    groups.setdefault(root, []).append(name)
                good_contigs_groups = [g for g in groups.values() if len(g) >= 2]
                for group in good_contigs_groups:
                    out.write(f"{sample}\t{gene}\t{','.join(group)}\n")
    # list of contigs to stitch for each gene

if __name__ == '__main__':
    main()