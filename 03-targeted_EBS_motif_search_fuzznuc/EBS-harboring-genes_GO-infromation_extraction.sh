#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: To identify the gene and its annotation for each promoter with EBS hit for the DEG set of candidate genes
# File requirement: 
# - fuzznuc outputs for the DEG sets
# - the modified genome .gff3 file (genes.gff3)
#///////////////////////////////////////////////////////////////////////
#
#
## 1@. ISOLATE THE LOCUS_TAG FOR EACH PROMOTER COORDINATES WITH EBS HIT
#############################################################################
# First, I remove the headings (there is a heading line for each gene promoter searched for EBSn motif)
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg1_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./05-Candidate-genes_GO-information/EBS-promG0hg1_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg2_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./05-Candidate-genes_GO-information/EBS-promG0hg2_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg3_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./05-Candidate-genes_GO-information/EBS-promG48hg3_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg4_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./05-Candidate-genes_GO-information/EBS-promG48hg4_noHeadings.txt
#
# Second, I isolate the first column with the coordinates of the hit within the promoter sequence, and I select only the start coordinate. I use uniq command to remove duplicates.
cat ./05-Candidate-genes_GO-information/EBS-promG0hg1_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./05-Candidate-genes_GO-information/EBS-promG0hg1_gene-coordinates.txt
cat ./05-Candidate-genes_GO-information/EBS-promG0hg2_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./05-Candidate-genes_GO-information/EBS-promG0hg2_gene-coordinates.txt
cat ./05-Candidate-genes_GO-information/EBS-promG48hg3_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./05-Candidate-genes_GO-information/EBS-promG48hg3_gene-coordinates.txt
cat ./05-Candidate-genes_GO-information/EBS-promG48hg4_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./05-Candidate-genes_GO-information/EBS-promG48hg4_gene-coordinates.txt
#
# Third, I use the start coordinate to extract the gene_name from the bed file. I use a for loop to do it for each coordinate...
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg1_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_G0hg1.bed >> ./05-Candidate-genes_GO-information/EBS-promG0hg1_geneName.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg2_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_G0hg2.bed >> ./05-Candidate-genes_GO-information/EBS-promG0hg2_geneName.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg3_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_G48hg3.bed >> ./05-Candidate-genes_GO-information/EBS-promG48hg3_geneName.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg4_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_G48hg4.bed >> ./05-Candidate-genes_GO-information/EBS-promG48hg4_geneName.txt; done
# ...and then I isolate it from the rest of columns and subsitute "MtrunA17chr" by "MtrunA17_chr" to create the gene_id.
cat ./05-Candidate-genes_GO-information/EBS-promG0hg1_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./05-Candidate-genes_GO-information/EBS-promG0hg1_geneID.txt
cat ./05-Candidate-genes_GO-information/EBS-promG0hg2_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./05-Candidate-genes_GO-information/EBS-promG0hg2_geneID.txt
cat ./05-Candidate-genes_GO-information/EBS-promG48hg3_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./05-Candidate-genes_GO-information/EBS-promG48hg3_geneID.txt
cat ./05-Candidate-genes_GO-information/EBS-promG48hg4_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./05-Candidate-genes_GO-information/EBS-promG48hg4_geneID.txt
#
## 2@. EXTRACT ANNOTATION INFORMATION FOR EACH GENE
#############################################################################
# Searching in the gene_name file for those genes that are already annotated in the genome using the locus_Tag file
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./05-Candidate-genes_GO-information/EBS-promG0hg1_annotation.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./05-Candidate-genes_GO-information/EBS-promG0hg2_annotation.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./05-Candidate-genes_GO-information/EBS-promG48hg3_annotation.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg4_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./05-Candidate-genes_GO-information/EBS-promG48hg4_annotation.txt; done
#
# Searching in the summary file where there is additional putative annotation information for non-characterized genes using the locus_Tag file
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./05-Candidate-genes_GO-information/EBS-promG0hg1_summary.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG0hg2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./05-Candidate-genes_GO-information/EBS-promG0hg2_summary.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./05-Candidate-genes_GO-information/EBS-promG48hg3_summary.txt; done
for i in $(cat ./05-Candidate-genes_GO-information/EBS-promG48hg4_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./05-Candidate-genes_GO-information/EBS-promG48hg4_summary.txt; done
#
## 3@. REMOVE UNNECESARY INTERMEDIATE FILES
#############################################################################
rm ./05-Candidate-genes_GO-information/*_noHeadings.txt 
rm ./05-Candidate-genes_GO-information/*_gene-coordinates.txt 
rm ./05-Candidate-genes_GO-information/*_geneName.txt 
rm ./05-Candidate-genes_GO-information/*_geneID.txt
#
# END