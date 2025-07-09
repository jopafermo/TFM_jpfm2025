#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: The goal of this script is to prepare the files required to perform the promoter isolation pipeline:
# - Gene file with only a column of gene_ids compatible with the locus_tag format of the gff3 file
# - Randomly selected non-DEG genes with similar sizes to the DEG files in the same 
# - Modified gff3 file with only the gene information (free of mRNA, CDS, and other features)
# Bioinformatic tools: bedtools (https://bedtools.readthedocs.io/en/latest/content/overview.html)
# File requirement: Any .tsv file harboring extracted genes from heat map clusters after edgeR analyses.
#///////////////////////////////////////////////////////////////////////
#
## 1@. FIRST: FORMATTING THE DEG TSV DATA SET TO OBTAIN A FILE WITH ONLY THE GENE_IDs
########################################################################################################################
# Create a file with the gene IDs of the upregulated genes in A17 at 0h and/or 48h from the heatmap analysis:
cat ./00-Selected_gene_list_files/DEG.down.G0h.III.hmtree.genes_1.tsv | sed 's/"/\t/g' | cut -f 2 > ./01-bed_files/G0h_hmG1.txt  # It has 27 genes
cat ./00-Selected_gene_list_files/DEG.down.G0h.III.hmtree.genes_2.tsv | sed 's/"/\t/g' | cut -f 2 > ./01-bed_files/G0h_hmG2.txt  # It has 8 genes
cat ./00-Selected_gene_list_files/DEG.down.G48h.III.hmtree.genes_3.tsv | sed 's/"/\t/g' | cut -f 2 > ./01-bed_files/G48h_hmG3.txt  # It has 172 genes
cat ./00-Selected_gene_list_files/DEG.down.G48h.III.hmtree.genes_4.tsv | sed 's/"/\t/g' | cut -f 2 > ./01-bed_files/G48h_hmG4.txt  # It has 141 genes
# this commands remove (""), select the column with the gene ids and remove the heading in the first row.
# I don't remove the (_) within the gene name, as it matches the locus_tag format of the gene names in the gff3 file.
#
## 2@. SECOND: FORMATTING THE NOSIG-DEG TSV DATA SET TO OBTAIN A FILE WITH ONLY THE GENE_IDs 
## & CREATING RANDOMLY SELECTED NON-DEG GENES SETS WITH SIMILAR SIZES TO THE DEG FILES
########################################################################################################################
# Create a file with the gene IDs of the differentially expressed genes with no statistical significance in M. truncatula roots 
# from the contrast skl_0h_vs_A17_0h.
cat ./00-Selected_gene_list_files/G0hI_NoSig_dataset.tsv | sed 's/"/\t/g' | cut -f4 > ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt  
#
# Create three replicates for four files of randomly selected NoSig-DEG with equivalent sizes to the DEG files
# Replicate 1
shuf -n 29 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_27DEG_rep1.txt
shuf -n 8 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_8DEG_rep1.txt
shuf -n 165 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_172DEG_rep1.txt
shuf -n 136 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_141DEG_rep1.txt
# Replicate 2
shuf -n 29 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_27DEG_rep2.txt
shuf -n 8 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_8DEG_rep2.txt
shuf -n 165 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_172DEG_rep2.txt
shuf -n 136 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_141DEG_rep2.txt
# Replicate 3
shuf -n 29 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_27DEG_rep3.txt
shuf -n 8 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_8DEG_rep3.txt
shuf -n 165 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_172DEG_rep3.txt
shuf -n 136 ./00-Selected_gene_list_files/nonSIG_DEG_ID.txt > ./00-Selected_gene_list_files/nonSig_141DEG_rep3.txt

## 3@. THIRD: I PREPARE A MODIFIED GFF3 FILE WITH ONLY THE GENE INFORMATION
#####################################################################
# In order to extract the matching lines for the gene ids obtained in the heatmap analysis,
# to have the coordinates of the genes, I want to get the gff3 file only with the gene information
# Since Column 2 (source) is not important for the analysis, and it contains the word "gene" I will remove it.
# And the empty columns 6 and 8 are not important for the analysis, but I will keep 6 for the bed 6 format score column
# and I will remove colunm 8 as I don't need it for the bed format
# The columns are: 1. seqid, 2. source, 3. type, 4. start, 5. end, 6. score, 7. strand, 8. phase, 9. attributes
# The command below will keep columns 1, 3, 4, 5, 7, and 9
# The output will be saved in a new file called modified_MtrunA17r5.0-ANR-EGN-r1.9.gff3
cat ./00-Selected_gene_list_files/MtrunA17r5.0-ANR-EGN-r1.9.gff3 | cut -f 1,3,4,5,6,7,9 > ./00-Selected_gene_list_files/modified_MtrunA17r5.0-ANR-EGN-r1.9.gff3
#
# Now I get the lines with "gene" pattern
# However, the lines with mRNA has the "gene" pattern too. So I use "ID=gene" for isolating only the gene lines.
grep "ID=gene" ././00-Selected_gene_list_files/modified_MtrunA17r5.0-ANR-EGN-r1.9.gff3 > ./01-bed_files/genes.gff3
# This file is the one I will use to extract the coordinates of the genes.
#
## END
