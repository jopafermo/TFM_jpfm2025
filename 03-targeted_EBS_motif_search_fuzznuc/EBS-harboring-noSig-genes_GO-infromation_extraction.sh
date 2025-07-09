#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: To identify the gene and its annotation for each promoter with EBS hit for the NonSIG DEG set of control genes
# File requirement: 
# - fuzznuc outputs for the NonSIG DEG sets
# - the modified genome .gff3 file (genes.gff3)
#///////////////////////////////////////////////////////////////////////
#
#
## 1@. ISOLATE THE LOCUS_TAG FOR EACH PROMOTER COORDINATES WITH EBS HIT
#############################################################################
# First, I remove the headings (there is a heading line for each gene promoter searched for EBSn motif)
# Repeat 1
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep1_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep1_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep1_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep1_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_noHeadings.txt
# Repeat 2
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep2_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep2_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep2_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep2_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_noHeadings.txt
# Repeat 3
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep3_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep3_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep3_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_noHeadings.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep3_filtered_EBSn.fuzznuc | grep -v "SeqName" > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_noHeadings.txt
#
# Second, I isolate the first column with the coordinates of the hit within the promoter sequence, and I select only the start coordinate. I use uniq command to remove duplicates.
# Repeat 1
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_gene-coordinates.txt
# Repeat 2
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_gene-coordinates.txt
# Repeat 3
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_gene-coordinates.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_noHeadings.txt | cut -f 1 | sed 's/-/\t/g' | cut -f 1 | uniq > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_gene-coordinates.txt
#
# Third, I use the start coordinate to extract the gene_name from the bed file. I use a for loop to do it for each coordinate...
# Repeat 1
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_8_rep1.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_27_rep1.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_172_rep1.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_141_rep1.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_geneName.txt; done
# Repeat 2
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_8_rep2.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_27_rep2.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_172_rep2.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_141_rep2.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_geneName.txt; done
# Repeat 3
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_8_rep3.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_27_rep3.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_172_rep3.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_geneName.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_gene-coordinates.txt); do grep -i $i ./02-promoters_files/2Kprom_noSigDEG_141_rep3.bed >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_geneName.txt; done
# ...and then I isolate it from the rest of columns and subsitute "MtrunA17chr" by "MtrunA17_chr" to create the gene_id.
# Repeat 1
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_geneID.txt
# Repeat 2
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_geneID.txt
# Repeat 3
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_geneID.txt
cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_geneName.txt | cut -f 4 | sed 's/MtrunA17/MtrunA17_/g' > ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_geneID.txt
#
## 2@. EXTRACT ANNOTATION INFORMATION FOR EACH GENE
#############################################################################
# Searching in the gene_name file for those genes that are already annotated in the genome using the locus_Tag file
# Repeat 1
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_annotation.txt; done
# Repeat 2
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_annotation.txt; done
# Repeat 3
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_annotation.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_annotation.txt; done
#
# Searching in the summary file where there is additional putative annotation information for non-characterized genes using the locus_Tag file
# Repeat 1
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep1_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep1_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep1_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep1_summary.txt; done
# Repeat 2
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep2_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep2_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep2_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep2_summary.txt; done
# Repeat 3
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_8_rep3_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_27_rep3_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_172_rep3_summary.txt; done
for i in $(cat ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_geneID.txt); do grep -i $i ./05-Candidate-genes_GO-information/MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv >> ./06-noSigDEG_withEBS_GOinformation/EBS-promnoSigDEG_141_rep3_summary.txt; done
#
## 3@. REMOVE UNNECESARY INTERMEDIATE FILES
#############################################################################
rm ./06-noSigDEG_withEBS_GOinformation/*_noHeadings.txt 
rm ./06-noSigDEG_withEBS_GOinformation/*_gene-coordinates.txt 
rm ./06-noSigDEG_withEBS_GOinformation/*_geneName.txt 
rm ./06-noSigDEG_withEBS_GOinformation/*_geneID.txt
#
# END