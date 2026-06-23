# ROOT <!-- omit in toc -->
- [Overview](#overview)
  - [How to use](#how-to-use)
  - [Dependencies](#dependencies)
  - [Citation](#citation)
- [1. Identifying Target Capture Homologous Gene Clusters within the Genomic Reference Dataset (Optional)](#1-identifying-target-capture-homologous-gene-clusters-within-the-genomic-reference-dataset-optional)
  - [1.1 Pull representative sequences from each homolog cluster](#11-pull-representative-sequences-from-each-homolog-cluster)
  - [1.2 Perform BLASTN](#12-perform-blastn)
  - [1.3 Sort BLASTN results](#13-sort-blastn-results)
- [2. Identify Paralogous Sequences for Multi-Copy HybPiper Assembly Reference](#2-identify-paralogous-sequences-for-multi-copy-hybpiper-assembly-reference)
  - [2.1 Make phylogenies](#21-make-phylogenies)
  - [2.2 Infer orthology and gene duplications](#22-infer-orthology-and-gene-duplications)
  - [2.3 Extract paralogous reference sequences](#23-extract-paralogous-reference-sequences)
- [3. Phylogenetically Group Orthologous Contigs](#3-phylogenetically-group-orthologous-contigs)
  - [3.1 Pull assembled hybpiper contigs and combine with genomic sequences](#31-pull-assembled-hybpiper-contigs-and-combine-with-genomic-sequences)
  - [3.2 Make phylogenies](#32-make-phylogenies)
  - [3.3 Group orthologous contigs from homolog gene trees](#33-group-orthologous-contigs-from-homolog-gene-trees)
  - [3.4 Stitch contigs and combine with reference sequences](#34-stitch-contigs-and-combine-with-reference-sequences)
  - [3.5 Make consensus species tree](#35-make-consensus-species-tree)

# Overview
This tree-based command-line workflow, **ROOT (Reference-guided Ortholog Organization of Target-capture contigs)**, prepares multi-copy reference to improve read mapping and separates HybPiper assembled target capture CDS contigs into orthologous groups, based on pre-existing genomic reference dataset. **The pre-existing reference dataset has to be genomic or transcriptomic** to infer gene duplications within your clade of interest. Target capture datasets cannot be used as the reference dataset in this workflow.

- **Step 1. Identifying Target Capture Homologous Gene Clusters within the Genomic Reference Dataset (Optional)**
  - This step is optional if you already have the corresponding homolog gene clusters from the reference dataset.
- **Step 2. Identify Paralogous Sequences for Multi-Copy HybPiper Assembly Reference**
  - This step produces a multi-copy reference that improves the read mapping and assembly of duplicated gene copies.
- **Step 3. Phylogenetically Group Orthologous Contigs**
  - This step phylogenetically groups HybPiper contigs into orthologous groups and stitches them together.

## How to use
This workflow is fully based on separate python scripts, no installation required, but many steps require outside software that needs installation (see dependencies). You can download the ROOT scripts by:

```console
git clone https://github.com/fkyeeb/ROOT.git
```
Or download the zip by clicking the green `Code` button on the top right.

## Dependencies
Some of the steps are optional, depending on the status of your pre-existing genomic reference dataset.
- BLAST (step 1.2) - [install instructions](https://www.ncbi.nlm.nih.gov/books/NBK569861/) or [with bioconda](https://anaconda.org/channels/bioconda/packages/blast/overview)
- parallel (step 2.1, 3.2, and 3.5, optional) - [install instructions](https://www.gnu.org/software/parallel/)
- mafft (step 2.1, 3.2, and 3.5) - [install instructions](https://mafft.cbrc.jp/alignment/software/source.html) or [with bioconda](https://anaconda.org/channels/bioconda/packages/mafft/overview)
- phyx (step 2.1, 3.2, and 3.5) - [install instructions](https://github.com/FePhyFoFum/phyx)
- raxml-ng (step 2.1, 3.2, and 3.5) - [install instructions](https://github.com/amkozlov/raxml-ng) or [with bioconda](https://anaconda.org/channels/bioconda/packages/raxml-ng/overview)
- H2O (step 2.2 and 3.5) - [install instructions](https://github.com/fkyeeb/H2O)
- ASTRAL (step 3.5) - [install instructions](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)

## Citation
If you used ROOT in your work, please cite:

- Feng K, Charboneau JLM, Smith SA. Beyond a Single Copy: Accurately Assemble Paralogs in Target Capture Datasets. In review.


Please also cite the publications of the outside softwares that you used in this workflow.

# 1. Identifying Target Capture Homologous Gene Clusters within the Genomic Reference Dataset (Optional)

This step identifies the homolog gene clusters in the genomic dataset that correspond to the target capture genes. If you already have the matching gene clusters, you can skip to step 2.
- **Input** - genomic homolog gene clusters (DNA sequences) in each individual fasta file in one designated folder
- **Output** - genomic homolog gene sequences corresponding to target capture genes, one fasta per gene

## 1.1 Pull representative sequences from each homolog cluster

```console
python pull_genomic_seqs.py <fasta_folder> <output_folder-optional> <sequence_#_per_cluster-optional>
```
- `<fasta_folder>` - where all the genomic homolog gene cluster fasta files are without other files.
- `<output_folder-optional>` - default is current folder if no argument is entered. There are two output files:
  - `genomic_seqs4blast.fa` - all the pulled genomic sequences.
  - `genomic_seq2cluster.tsv` - lists all pulled genomic sequences to its corresponding cluster name.
- `<sequence_#_per_cluster-optional>` - default is 2, if no argument is entered. It has to be an integer.

## 1.2 Perform BLASTN
First make the blast database with `makeblastdb` and perform `blastn` to report all hits with evalue < 1e-6.

```console
makeblastdb -in genomic_seqs4blast.fa -parse_seqids -dbtype nucl -out <preferred_id>
blastn -query <target_capture_reference_genes> -db <folder>/<preferred_id> -task dc-megablast -evalue 1e-6 -outfmt "6 qseqid sseqid pident qlen slen sstart send length mismatch gapopen evalue" -out <output_folder>/<output_file_name>
```
- `<target_capture_reference_genes>` - one ingroup sequence per gene. If there is no ingroup, then the closest outgroup sequence. Reference sequence name has to be `sample-gene_name`. `gene_name` should not contain "-".
- `-outfmt` - this output format is required to process the blast result. If you used a different format, you would need to manually change the blast result processing script `sort_blast_result.py` to fit your result format.

## 1.3 Sort BLASTN results
This script filters through the blast results and ignores matches when the alignment with the query (reference) sequence is <90% of its full length **AND** the alignment with the target (genomic) sequence is <50% of its full length. Then, write the matched genomic homolog clusters to new fasta files in the corresponding target capture gene names.

```console
python sort_blast_result.py <blast_result_file> <genomic_seq2cluster.tsv_location> <all_or_best> <fasta_input_folder> <fasta_output_folder-optional> <tsv_output_folder-optional>
```
- `<blast_result_file>` and `<genomic_seq2cluster.tsv_location>` - are the folder location + the file names.
- `<all_or_best>` - If "all", all matched homolog clusters will be written into the output fasta file. If "best", only the best hit cluster, i.e., with the longest alignment with the query, will be written into the output fasta file. If you want to be conservative, `best` is recommended. `all` can increase taxon sampling but can potentially include non-homologous clusters.
- `<fasta_input_folder>` - where all the genomic homolog gene cluster fasta files are without other files.
- `<fasta_output_folder-optional>` - default is current folder if no argument is entered. It is where the matched genomic homolog clusters will be written to, in the corresponding target capture gene names.
- `<tsv_output_folder-optional>` - default is current folder if no argument is entered. `blast_result_filtered.tsv` reports the `gene cluster genomic_seq qry_perc sub_perc`. "qry_perc" is the percentage of alignment of the query sequence; "sub_perc" is that of the target sequence.

# 2. Identify Paralogous Sequences for Multi-Copy HybPiper Assembly Reference

As ROOT manuscript documented, including sequences from all inferred duplicated clades in HybPiper reference may improve the recovery of paralogous sequences. Thus, here we present a way to retrieve paralogous sequences from the genomic reference dataset for HybPiper reference. 

- **Input** - genomic homolog gene fasta corresponding to target capture genes in one designated folder
- **Output** - reference sequences fasta file for HybPiper assembly

## 2.1 Make phylogenies
You don't need to follow exactly these steps if you already have a confortable pipeline to make trees. You can also skip this step if you already have the genomic homolog gene trees ready.

Some of these steps might take a few hours if you do not parallelize the process, so we recommend using `parallel` but you don't have to.

**2.1.1 Align the fasta files**

Feel free to use a different alignment tool and change the commands.

- **Input** - genomic homolog gene sequence fasta corresponding to target capture genes. The gene sequence name format needs to be `sample_id@seq_id` to prepare for the future steps. 
- **Output** - aligned fasta (*.aln)

With `parallel`
```console
parallel -j <#_of_jobs> "mafft --auto --maxiterate 1000 --thread <#_of_threads> {} > {.}.aln" ::: *.<fasta_file_ending>
```
- `<#_of_jobs>` * `<#_of_threads>` - cannot be more than the number of threads of your machine. It should be slightly less than that.
- `<fasta_file_ending>` - e.g. "fa"

Without `parallel`
```console
for f in *.<fasta_file_ending>; do
    mafft --auto --maxiterate 1000 --thread <#_of_threads> "$f" > "${f%.fa}.aln"
done
```

**2.1.2 Remove sites that are mostly gaps across sample**

- **Input** - aligned fasta (*.aln)
- **Output** - cleaned alignment (*.aln-cln)

With `parallel`
```console
parallel -j <#_of_jobs> "pxclsq -s {} -p <propn_data_present> -o {}-cln" ::: *.aln
```
- `<#_of_jobs>` - cannot be more than the number of threads of your machine. It should be slightly less than that.
- `<propn_data_present>` - for each site, the minimum proportion of samples that have a base present. Site with a smaller proportion will be removed. We used `0.3` in the ROOT manuscript. Up to your judgment!

Without `parallel`
```console
for f in *.aln; do
    pxclsq -s "$f" -p <propn_data_present> -o "${f}-cln"
done
```

**2.1.3 Make phylogenies**

Feel free to use a different phylogenetic tool and change the commands.

- **Input** - cleaned alignment (*.aln-cln)
- **Output** - phylogenetic trees (*.bestTree)

With `parallel`
```console
parallel -j <#_of_jobs> "raxml-ng --msa {} --threads <#_of_threads> --model GTR+G" ::: *.aln-cln
```
- `<#_of_jobs>` * `<#_of_threads>` - cannot be more than the number of threads of your machine. It should be slightly less than that.

Without `parallel`
```console
for f in *.aln-cln; do
    raxml-ng --msa "$f" --threads <#_of_threads> --model GTR+G
done
```

## 2.2 Infer orthology and gene duplications

We need to infer gene orthology/paralogy before we can come up with the list of references. More details of `h2o` commands see [tutorial](https://github.com/fkyeeb/H2O/blob/main/tutorials/tutorial.md) of the tool. Note that `h2o` will skip trees when the outgroups are polyphyletic or nonexistent. 

- **Input** - homolog gene trees
- **Output** - processed gene trees with gene duplications identified as node labels (processed_trees/unpruned/*_rooted_processed.tre)

```console
h2o infer_ortho -t <homolog_tree_folder> -e <tree_file_ending> -f <id2species_name_file> -of <outgroup_file> > infer_ortho.log
```
- `-f <id2species_name_file>` - your tip names should be unique sequence ID formated as `sample_id@seq_id`. And this `sample_id@seq_id` name should match exactly with the names in the homolog gene fasta files. Each line in the file should be `sample_id<tab>preferred_tip_name`. 
- `<tree_file_ending>` - e.g. ".tre", ".bestTree"

> [!TIP]
> If you have multiple samples per some species and you do not want to account for intraspecific duplications, please give them the same preferred tip name in `<id2species_name_file>`. `h2o` will take them as separate taxa if you do not give them the same name.

## 2.3 Extract paralogous reference sequences

Make sure `tree_reader.py` and `node.py` are in the same folder that you are running this command. Each gene will have at least 3 (default, or `<#_of_ref_per_ortho-optional>`) longest sequences from the cluster as reference. If there are inferred gene duplications, 1~3 (default, or `<#_of_ref_per_ortho-optional>`) paralog sequences will be added to the reference for each duplicated clade. The number depends on the size of the duplicated clade. If the tree was skipped by `h2o`, 3 (default, or `<#_of_ref_per_ortho-optional>`) longest sequences will be selected as references.

```console
python select_paralog_reference.py <outgroup_file> <h2o_processed_trees_folder> <fasta_folder> <gene_name_file> <fasta_file_ending> <#_of_ref_per_ortho-optional> <reference_file_name-optional>
```
- `<outgroup_file>` - one outgroup sample id (in `sample_id@seq_id`) per line. They can be outgroups or whichever samples you do not want to include in the references.
- `<h2o_processed_trees_folder>` - is the location of the `processed_trees/` folder. Enter `processed_trees/` if the folder is in the currect directry.
- `<fasta_folder>` - folder containing homolog gene sequence fasta files, file names formatted as `<gene_name>.<fasta_file_ending>`
- `<gene_names_file>` - all gene names in a file, one name per line.
- `<fasta_file_ending>` - e.g. `fa`
- `<#_of_ref_per_ortho-optional>` - the number of sequences to use as reference per ortho group. Default is 3 if no entry.
- `<reference_file_name-optional>` - default is `hybpiper_paralog_reference.fa` if no entry.

> [!TIP]
> In the h2o `processed_trees/unpruned/` folder, the processed homolog tree names should be `<gene_name>.*rooted_processed_id.tre` and the fasta file names in the `<fasta_folder>` should be `<gene_name>.<fasta_file_ending>`. If not, the script will not be able to open the the correct corresponding fasta file.

# 3. Phylogenetically Group Orthologous Contigs

- **Input** - hybpiper assembled contigs
- **Output** - orthologous stitch contigs (to step 3.4), or to consensus tree (to step 3.5)

## 3.1 Pull assembled hybpiper contigs and combine with genomic sequences

Output fasta files of hybpiper contigs and genomic sequences for each gene in each sample. Each sample has its own folder in the output folder. Output fasta files are written into each sample's folder. Also output `all_contigs.fa` and `all_contig_positions.tsv` in the output directory. `all_contigs.fa` contains all the filtered contigs and `all_contig_positions.tsv` records all alignment information of the contigs with the reference sequences, which will be used for contig stitching. If the gene only has one contig for this sample, the fasta file will not be written to the sample folder, because it doesn't need to be grouped. But this contig will go into `all_contigs.fa`.

```console
python extract_hybpiper_contigs.py <sample_names_file> <gene_names_file> <genomic_sequence_folder> <hybpiper_result_folder> <output_folder> <fasta_file_ending-optional>
```
- `<sample_names_file>` - all sample names in a file, one name per line.
- `<gene_names_file>` - all gene names in a file, one name per line.
- `<genomic_sequence_folder>` - folder for genomic sequence corresponding to each gene, file names should be `<gene_name>.fa`, unless other `fasta_file_ending` specified.
- `<hybpiper_result_folder>` - the hybpiper base directory containing the hybpiper result files and gene folders. Enter `.` if it is the current directory.
- `<fasta_file_ending-optional>` - default is `fa`, if no entry. Do not include `.`, example entry: `fasta`.

## 3.2 Make phylogenies

Repeat [2.1 Make phylogenies](#21-make-phylogenies) with the output fasta files from step 3.1. You can use alternative methods if you have a different phylogenetic pipeline that you are comfortable with. The only thing is different here is adding constraint trees to raxml, because adding the contigs can potentially alter the relationships among reference sequences. Adding the constraint tree can also speed up the phylogenetic reconstructions. The altered commands are:

With `parallel`
```console
parallel -j <#_of_jobs> "raxml-ng --msa {} --threads <#_of_threads> --model GTR+G --tree-constraint <constraint_tree>" ::: *.aln-cln
```
- `<#_of_jobs>` * `<#_of_threads>` - cannot be more than the number of threads of your machine. It should be slightly less than that.
- `<constraint_tree>` should be entered as `/path/to/file/{.}.tree_file_ending`, e.g. `ROOT/{.}.aln-cln.raxml.bestTree`. The constraint trees should be the ones constructed from step 2.1.3

Without `parallel`
```console
for f in *.aln-cln; do
    raxml-ng --msa "$f" --threads <#_of_threads> --model GTR+G --tree-constraint <constraint_tree>
done
```

## 3.3 Group orthologous contigs from homolog gene trees

Output a `contig_groups-<tree_file_folder>.tsv` file contains all ortholog contig group information in the current directory. If this gene has no valid group, it will not be listed in this file. This script will also print a percentage of total contigs grouped (for genes with more than one contig), which may help verify the choice of NN and BL thresholds.

```console
python group_contigs.py <tree_file_folder> <tree_file_ending> <sample_names_file> <NN-optional> <BL-optional>
```
- `<tree_file_folder>` - the folder containing the phylogenies from step 3.2
- `<tree_file_ending>` - `.aln-cln.raxml.bestTree` if you followed this tutorial. The script expects the tree file names to be `<gene_name><tree_file_ending>`.
- `<sample_names_file>` - all sample names in a file, one name per line.
- `<NN-optional>` - the number of nodes threshold for grouping contigs. Default is 5 if no argument is entered. The cutoff should be applicable to other datasets, but you should check the reported percentage and trees for sanity check.
- `<BL-optional>` - the branch length threshold for grouping contigs. Default is 0.05 if no argument is entered. The cutoff should be applicable to other datasets, but you should check the reported percentage and trees for sanity check.

## 3.4 Stitch contigs and combine with reference sequences

First, stitches the contigs together. In the output folder, the stitched contigs of each gene is written in `*_stitched.fa`. Also, `stitched_all.fa` contains all stitched contigs, `stitched_info.tsv` records all the contigs are stitched together, and `copy_counts.tsv` records how many copies of each gene is produced for each sample, written in the output folder. If a gene has no valid group but has contigs, it will output the longest contig.
```console
python stitch_contigs.py <sample_names_file> <gene_names_file> <contig_folder> <contig_groups_tsv> <output_folder> <max_copies-optional> <min_overlap-optional>
```
- `<contig_folder>` - is the `<output_folder>` in step 3.1
- `<contig_groups_tsv>` -  is the output of step 3.2
- `<max_copies-optional>` -  is the max gene copy to output for each gene. Default is no limit if no argument is entered.
- `<min_overlap-optional>` - is the min bp overlap between stitched contigs.  Default is 50 if no argument is entered. If one gene has more than one group of contigs, and the stitched contigs from them do not have >50 bp overlap in alignment, the smaller groups will be dropped. Because the grouping is potentially due to intragenic conflict instead of paralogy.

Then, combine the stitched contigs with the reference sequence, preparing for phylogenetic reconstruction. Combined fasta files will be written into the output folder. This step is optional if you do not want to make individual gene trees.
```console
python combine_seqs.py <sample_names_file> <gene_names_file> <ref_fasta_folder> <stitched_contig_folder> <output_folder> <ref_fasta_ending-optional>
```
-  `<ref_fasta_ending-optional>` - Default is `fa` if no argument is entered. The full fasta file names should be in the format of `<gene>.<ref_fasta_ending>`

## 3.5 Make consensus species tree

If you would like to make a consensus tree with ASTRAL, you can repeat step 3.2 and 2.1 to make homologous gene trees and infer orthology.Then, in the `h2o` result folder, make input file for ASTRAL, run it, and root the output tree.
```console
cat $(ls processed_trees/unpruned/*ortho*.tre | grep -v '_id\.tre') > ASTRAL_in.tre
astral4 -i ASTRAL_in.tre -o ASTRAL_out.tre -t <#_of_threads>
pxrr -t ASTRAL_out.tre -f <outgroup_file> -o ASTRAL_rooted.tre
```
- `<outgroup_file>` - one outgroup sample name per line.