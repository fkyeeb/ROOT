import os,sys

args = sys.argv
if len(args) < 2:
    print("Usage: python pull_genomic_seqs.py <fasta_folder> <output_folder-optional> <sequence_#_per_cluster-optional>")
    sys.exit(1)
folder = args[1]

if len(args) > 3:
    num = int(args[3])
else:
    num = 2

if len(args) > 2:
    if args[2][-1] != "/":
        args[2] += "/"
    output = args[2] + "genomic_seqs4blast.fa"
    output_sum = args[2] + "genomic_seq2cluster.tsv"
else:
    output = "genomic_seqs4blast.fa"
    output_sum = "genomic_seq2cluster.tsv"


output_f = open(output,"w")
output_sum_f = open(output_sum,"w")
count = 0

for file in os.listdir(folder):
    if file.endswith(".fa"):
        path = os.path.join(folder, file)
        # cluster = file.split("_")[0]
        cluster = file
        with open(path,"r") as f:
            dic = {}
            seq = ""
            for line in f:
                if line.startswith(">"):
                    if seq != "":
                        dic[name] = seq + "\n"
                        seq = ""
                    name = line
                else:
                    seq += line.strip()
            dic[name] = seq + "\n"
    dic = dict(sorted(dic.items(), key=lambda x: len(x[1]), reverse=True))
    for key, value in list(dic.items())[:num]:
        output_sum_f.write(key.strip()[1:] + "\t" + cluster + "\n")
        output_f.write(key)
        output_f.write(value)
    count += 1

output_f.close()
print(f"Finished processing {count} clusters.")