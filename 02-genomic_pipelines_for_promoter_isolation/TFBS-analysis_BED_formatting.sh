#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: This script is used to obtain the BED6 format files for both the DEG genes and the nonSIG-DEG genes.
# Bioinformatic tools: bedtools (https://bedtools.readthedocs.io/en/latest/content/overview.html)
# File required: this script will uses the output from the previous script
# - TFBS-analysis_DATA_prep.sh
#///////////////////////////////////////////////////////////////////////
#
## 1@. FIRST: EXTRACTING PROMOTORES OF THE DE GENE SETS AND CREATING THE BED6 FILES
##########################################################################################
# Extracting the lines matching the locus_tags.
# With -A0 I get only the matching lines
grep -A0 -f ./01-bed_files/G0h_hmG1.txt ./01-bed_files/genes.gff3 > ./01-bed_files/G0h_hmG1_results.txt 
grep -A0 -f ./01-bed_files/G0h_hmG2.txt ./01-bed_files/genes.gff3 > ./01-bed_files/G0h_hmG2_results.txt
grep -A0 -f ./01-bed_files/G48h_hmG3.txt ./01-bed_files/genes.gff3 > ./01-bed_files/G48h_hmG3_results.txt  
grep -A0 -f ./01-bed_files/G48h_hmG4.txt ./01-bed_files/genes.gff3 > ./01-bed_files/G48h_hmG4_results.txt
#
# In order to create the BED6 file, I have to rearrange the columns. 
# So, I created individual files for each column and I pasted them together in the right order.
# @1. Column with chromosome names
cut -f1 ./01-bed_files/G0h_hmG1_results.txt > ./01-bed_files/chrom_G0hg1.txt
cut -f1 ./01-bed_files/G0h_hmG2_results.txt > ./01-bed_files/chrom_G0hg2.txt
cut -f1 ./01-bed_files/G48h_hmG3_results.txt > ./01-bed_files/chrom_G48hg3.txt
cut -f1 ./01-bed_files/G48h_hmG4_results.txt > ./01-bed_files/chrom_G48hg4.txt
# @2. Columns with start and end coordinates
cut -f3,4 ./01-bed_files/G0h_hmG1_results.txt  > ./01-bed_files/coordinates_G0hg1.txt
cut -f3,4 ./01-bed_files/G0h_hmG2_results.txt  > ./01-bed_files/coordinates_G0hg2.txt
cut -f3,4 ./01-bed_files/G48h_hmG3_results.txt  > ./01-bed_files/coordinates_G48hg3.txt
cut -f3,4 ./01-bed_files/G48h_hmG4_results.txt  > ./01-bed_files/coordinates_G48hg4.txt
# @3. Column with score information (can be empty)
cut -f5 ./01-bed_files/G0h_hmG1_results.txt > ./01-bed_files/empty_score_G0hg1.txt 
cut -f5 ./01-bed_files/G0h_hmG2_results.txt > ./01-bed_files/empty_score_G0hg2.txt
cut -f5 ./01-bed_files/G48h_hmG3_results.txt > ./01-bed_files/empty_score_G48hg3.txt 
cut -f5 ./01-bed_files/G48h_hmG4_results.txt > ./01-bed_files/empty_score_G48hg4.txt
# @4. Column with strand information
cut -f6 ./01-bed_files/G0h_hmG1_results.txt > ./01-bed_files/strand_G0hg1.txt
cut -f6 ./01-bed_files/G0h_hmG2_results.txt > ./01-bed_files/strand_G0hg2.txt
cut -f6 ./01-bed_files/G48h_hmG3_results.txt  > ./01-bed_files/strand_G48hg3.txt
cut -f6 ./01-bed_files/G48h_hmG4_results.txt  > ./01-bed_files/strand_G48hg4.txt
# @5. Column with gene IDs and other attributes
cut -f7 ./01-bed_files/G0h_hmG1_results.txt > ./01-bed_files/IDinfo_G0hg1.txt
cut -f7 ./01-bed_files/G0h_hmG2_results.txt > ./01-bed_files/IDinfo_G0hg2.txt 
cut -f7 ./01-bed_files/G48h_hmG3_results.txt > ./01-bed_files/IDinfo_G48hg3.txt
cut -f7 ./01-bed_files/G48h_hmG4_results.txt > ./01-bed_files/IDinfo_G48hg4.txt 
# @6. Isolating only the gene names as a single column
sed 's/;/\t/g' < ./01-bed_files/IDinfo_G0hg1.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_G0hg1.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo_G0hg2.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_G0hg2.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_G48hg3.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_G48hg3.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo_G48hg4.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_G48hg4.txt 
# @7. Paste the columns together to create the BED6 file
paste ./01-bed_files/chrom_G0hg1.txt ./01-bed_files/coordinates_G0hg1.txt ./01-bed_files/geneNames_G0hg1.txt ./01-bed_files/empty_score_G0hg1.txt ./01-bed_files/strand_G0hg1.txt > ./01-bed_files/A17GeneUpG0hg1.bed
paste ./01-bed_files/chrom_G0hg2.txt ./01-bed_files/coordinates_G0hg2.txt ./01-bed_files/geneNames_G0hg2.txt ./01-bed_files/empty_score_G0hg2.txt ./01-bed_files/strand_G0hg2.txt > ./01-bed_files/A17GeneUpG0hg2.bed
paste ./01-bed_files/chrom_G48hg3.txt ./01-bed_files/coordinates_G48hg3.txt ./01-bed_files/geneNames_G48hg3.txt ./01-bed_files/empty_score_G48hg3.txt ./01-bed_files/strand_G48hg3.txt > ./01-bed_files/A17GeneUpG48hg3.bed
paste ./01-bed_files/chrom_G48hg4.txt ./01-bed_files/coordinates_G48hg4.txt ./01-bed_files/geneNames_G48hg4.txt ./01-bed_files/empty_score_G48hg4.txt ./01-bed_files/strand_G48hg4.txt > ./01-bed_files/A17GeneUpG48hg4.bed
# @8. Removing the intermediate files
rm ./01-bed_files/G0h_hmG1_results.txt ./01-bed_files/G0h_hmG2_results.txt ./01-bed_files/G48h_hmG3_results.txt ./01-bed_files/G48h_hmG4_results.txt
rm ./01-bed_files/chrom_G0hg1.txt ./01-bed_files/chrom_G0hg2.txt ./01-bed_files/chrom_G48hg3.txt ./01-bed_files/chrom_G48hg4.txt
rm ./01-bed_files/coordinates_G0hg1.txt ./01-bed_files/coordinates_G0hg2.txt ./01-bed_files/coordinates_G48hg3.txt ./01-bed_files/coordinates_G48hg4.txt
rm ./01-bed_files/empty_score_G0hg1.txt ./01-bed_files/empty_score_G0hg2.txt ./01-bed_files/empty_score_G48hg3.txt ./01-bed_files/empty_score_G48hg4.txt
rm ./01-bed_files/strand_G0hg1.txt ./01-bed_files/strand_G0hg2.txt ./01-bed_files/strand_G48hg3.txt ./01-bed_files/strand_G48hg4.txt
rm ./01-bed_files/IDinfo_G0hg1.txt ./01-bed_files/IDinfo_G0hg2.txt ./01-bed_files/IDinfo_G48hg3.txt ./01-bed_files/IDinfo_G48hg4.txt
rm ./01-bed_files/geneNames_G0hg1.txt ./01-bed_files/geneNames_G0hg2.txt ./01-bed_files/geneNames_G48hg3.txt ./01-bed_files/geneNames_G48hg4.txt
#
## 1@. FIRST: EXTRACTING PROMOTORES OF THE DE GENE SETS AND CREATING THE BED6 FILES
##########################################################################################
# Extracting the lines matching the locus_tags.
# With -A0 I get only the matching lines
# Replicate 1
grep -A0 -f ./00-Selected_gene_list_files/nonSig_27DEG_rep1.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_27DEG_rep1_results.txt  
grep -A0 -f ./00-Selected_gene_list_files/nonSig_8DEG_rep1.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_8DEG_rep1_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_172DEG_rep1.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_172DEG_rep1_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_141DEG_rep1.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_141DEG_rep1_results.txt
# Replicate 2
grep -A0 -f ./00-Selected_gene_list_files/nonSig_27DEG_rep2.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_27DEG_rep2_results.txt  
grep -A0 -f ./00-Selected_gene_list_files/nonSig_8DEG_rep2.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_8DEG_rep2_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_172DEG_rep2.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_172DEG_rep2_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_141DEG_rep2.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_141DEG_rep2_results.txt
# Replicate 3
grep -A0 -f ./00-Selected_gene_list_files/nonSig_27DEG_rep3.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_27DEG_rep3_results.txt  
grep -A0 -f ./00-Selected_gene_list_files/nonSig_8DEG_rep3.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_8DEG_rep3_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_172DEG_rep3.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_172DEG_rep3_results.txt
grep -A0 -f ./00-Selected_gene_list_files/nonSig_141DEG_rep3.txt ./01-bed_files/genes.gff3 > ./01-bed_files/nonSig_141DEG_rep3_results.txt
#
# In order to create the BED6 file, I have to rearrange the columns. 
# So, I created individual files for each column and I pasted them together in the right order.
# @1. Column with chromosome names
# Replicate 1
cut -f1 ./01-bed_files/nonSig_27DEG_rep1_results.txt > ./01-bed_files/chrom_1.txt
cut -f1 ./01-bed_files/nonSig_8DEG_rep1_results.txt > ./01-bed_files/chrom_2.txt
cut -f1 ./01-bed_files/nonSig_172DEG_rep1_results.txt > ./01-bed_files/chrom_3.txt
cut -f1 ./01-bed_files/nonSig_141DEG_rep1_results.txt > ./01-bed_files/chrom_4.txt
# Replicate 1
cut -f1 ./01-bed_files/nonSig_27DEG_rep2_results.txt > ./01-bed_files/chrom_5.txt
cut -f1 ./01-bed_files/nonSig_8DEG_rep2_results.txt > ./01-bed_files/chrom_6.txt
cut -f1 ./01-bed_files/nonSig_172DEG_rep2_results.txt > ./01-bed_files/chrom_7.txt
cut -f1 ./01-bed_files/nonSig_141DEG_rep2_results.txt > ./01-bed_files/chrom_8.txt
# Replicate 1
cut -f1 ./01-bed_files/nonSig_27DEG_rep3_results.txt > ./01-bed_files/chrom_9.txt
cut -f1 ./01-bed_files/nonSig_8DEG_rep3_results.txt > ./01-bed_files/chrom_10.txt
cut -f1 ./01-bed_files/nonSig_172DEG_rep3_results.txt > ./01-bed_files/chrom_11.txt
cut -f1 ./01-bed_files/nonSig_141DEG_rep3_results.txt > ./01-bed_files/chrom_12.txt
# @2. Columns with start and end coordinates
# Replicate 1
cut -f3,4 ./01-bed_files/nonSig_27DEG_rep1_results.txt  > ./01-bed_files/coordinates_1.txt
cut -f3,4 ./01-bed_files/nonSig_8DEG_rep1_results.txt  > ./01-bed_files/coordinates_2.txt
cut -f3,4 ./01-bed_files/nonSig_172DEG_rep1_results.txt  > ./01-bed_files/coordinates_3.txt
cut -f3,4 ./01-bed_files/nonSig_141DEG_rep1_results.txt  > ./01-bed_files/coordinates_4.txt
# Replicate 2
cut -f3,4 ./01-bed_files/nonSig_27DEG_rep2_results.txt  > ./01-bed_files/coordinates_5.txt
cut -f3,4 ./01-bed_files/nonSig_8DEG_rep2_results.txt  > ./01-bed_files/coordinates_6.txt
cut -f3,4 ./01-bed_files/nonSig_172DEG_rep2_results.txt  > ./01-bed_files/coordinates_7.txt
cut -f3,4 ./01-bed_files/nonSig_141DEG_rep2_results.txt  > ./01-bed_files/coordinates_8.txt
# Replicate 3
cut -f3,4 ./01-bed_files/nonSig_27DEG_rep3_results.txt  > ./01-bed_files/coordinates_9.txt
cut -f3,4 ./01-bed_files/nonSig_8DEG_rep3_results.txt  > ./01-bed_files/coordinates_10.txt
cut -f3,4 ./01-bed_files/nonSig_172DEG_rep3_results.txt  > ./01-bed_files/coordinates_11.txt
cut -f3,4 ./01-bed_files/nonSig_141DEG_rep3_results.txt  > ./01-bed_files/coordinates_12.txt
# @3. Column with score information (can be empty)
# Replicate 1
cut -f5 ./01-bed_files/nonSig_27DEG_rep1_results.txt > ./01-bed_files/empty_score_1.txt 
cut -f5 ./01-bed_files/nonSig_8DEG_rep1_results.txt > ./01-bed_files/empty_score_2.txt
cut -f5 ./01-bed_files/nonSig_172DEG_rep1_results.txt > ./01-bed_files/empty_score_3.txt
cut -f5 ./01-bed_files/nonSig_141DEG_rep1_results.txt > ./01-bed_files/empty_score_4.txt
# Replicate 2
cut -f5 ./01-bed_files/nonSig_27DEG_rep2_results.txt > ./01-bed_files/empty_score_5.txt 
cut -f5 ./01-bed_files/nonSig_8DEG_rep2_results.txt > ./01-bed_files/empty_score_6.txt
cut -f5 ./01-bed_files/nonSig_172DEG_rep2_results.txt > ./01-bed_files/empty_score_7.txt
cut -f5 ./01-bed_files/nonSig_141DEG_rep2_results.txt > ./01-bed_files/empty_score_8.txt
# Replicate 3
cut -f5 ./01-bed_files/nonSig_27DEG_rep3_results.txt > ./01-bed_files/empty_score_9.txt 
cut -f5 ./01-bed_files/nonSig_8DEG_rep3_results.txt > ./01-bed_files/empty_score_10.txt
cut -f5 ./01-bed_files/nonSig_172DEG_rep3_results.txt > ./01-bed_files/empty_score_11.txt
cut -f5 ./01-bed_files/nonSig_141DEG_rep3_results.txt > ./01-bed_files/empty_score_12.txt
# @4. Column with strand information
# Replicate 1
cut -f6 ./01-bed_files/nonSig_27DEG_rep1_results.txt > ./01-bed_files/strand_1.txt
cut -f6 ./01-bed_files/nonSig_8DEG_rep1_results.txt > ./01-bed_files/strand_2.txt
cut -f6 ./01-bed_files/nonSig_172DEG_rep1_results.txt > ./01-bed_files/strand_3.txt
cut -f6 ./01-bed_files/nonSig_141DEG_rep1_results.txt > ./01-bed_files/strand_4.txt
# Replicate 2
cut -f6 ./01-bed_files/nonSig_27DEG_rep2_results.txt > ./01-bed_files/strand_5.txt
cut -f6 ./01-bed_files/nonSig_8DEG_rep2_results.txt > ./01-bed_files/strand_6.txt
cut -f6 ./01-bed_files/nonSig_172DEG_rep2_results.txt > ./01-bed_files/strand_7.txt
cut -f6 ./01-bed_files/nonSig_141DEG_rep2_results.txt > ./01-bed_files/strand_8.txt
# Replicate 3
cut -f6 ./01-bed_files/nonSig_27DEG_rep3_results.txt > ./01-bed_files/strand_9.txt
cut -f6 ./01-bed_files/nonSig_8DEG_rep3_results.txt > ./01-bed_files/strand_10.txt
cut -f6 ./01-bed_files/nonSig_172DEG_rep3_results.txt > ./01-bed_files/strand_11.txt
cut -f6 ./01-bed_files/nonSig_141DEG_rep3_results.txt > ./01-bed_files/strand_12.txt
# @5. Column with gene IDs and other attributes
# Replicate 1
cut -f7 ./01-bed_files/nonSig_27DEG_rep1_results.txt > ./01-bed_files/IDinfo_1.txt
cut -f7 ./01-bed_files/nonSig_8DEG_rep1_results.txt > ./01-bed_files/IDinfo_2.txt 
cut -f7 ./01-bed_files/nonSig_172DEG_rep1_results.txt > ./01-bed_files/IDinfo_3.txt
cut -f7 ./01-bed_files/nonSig_141DEG_rep1_results.txt > ./01-bed_files/IDinfo_4.txt 
# Replicate 2
cut -f7 ./01-bed_files/nonSig_27DEG_rep2_results.txt > ./01-bed_files/IDinfo_5.txt
cut -f7 ./01-bed_files/nonSig_8DEG_rep2_results.txt > ./01-bed_files/IDinfo_6.txt 
cut -f7 ./01-bed_files/nonSig_172DEG_rep2_results.txt > ./01-bed_files/IDinfo_7.txt
cut -f7 ./01-bed_files/nonSig_141DEG_rep2_results.txt > ./01-bed_files/IDinfo_8.txt 
# Replicate 3
cut -f7 ./01-bed_files/nonSig_27DEG_rep3_results.txt > ./01-bed_files/IDinfo_9.txt
cut -f7 ./01-bed_files/nonSig_8DEG_rep3_results.txt > ./01-bed_files/IDinfo_10.txt 
cut -f7 ./01-bed_files/nonSig_172DEG_rep3_results.txt > ./01-bed_files/IDinfo_11.txt
cut -f7 ./01-bed_files/nonSig_141DEG_rep3_results.txt > ./01-bed_files/IDinfo_12.txt 
# @6. Isolating only the gene names as a single column
# Replicate 1
sed 's/;/\t/g' < ./01-bed_files/IDinfo_1.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_1.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo_2.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_2.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_3.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_3.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_4.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_4.txt 
# Replicate 2
sed 's/;/\t/g' < ./01-bed_files/IDinfo_5.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_5.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo_6.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_6.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_7.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_7.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_8.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_8.txt 
# Replicate 3
sed 's/;/\t/g' < ./01-bed_files/IDinfo_9.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_9.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo_10.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_10.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_11.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_11.txt 
sed 's/;/\t/g' < ./01-bed_files/IDinfo_12.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames_12.txt 
# @7. Paste the columns together to create the BED6 file
# Replicate 1
paste ./01-bed_files/chrom_1.txt ./01-bed_files/coordinates_1.txt ./01-bed_files/geneNames_1.txt ./01-bed_files/empty_score_1.txt ./01-bed_files/strand_1.txt > ./01-bed_files/nonSigGenes_27_rep1.bed
paste ./01-bed_files/chrom_2.txt ./01-bed_files/coordinates_2.txt ./01-bed_files/geneNames_2.txt ./01-bed_files/empty_score_2.txt ./01-bed_files/strand_2.txt > ./01-bed_files/nonSigGenes_8_rep1.bed
paste ./01-bed_files/chrom_3.txt ./01-bed_files/coordinates_3.txt ./01-bed_files/geneNames_3.txt ./01-bed_files/empty_score_3.txt ./01-bed_files/strand_3.txt > ./01-bed_files/nonSigGenes_172_rep1.bed
paste ./01-bed_files/chrom_4.txt ./01-bed_files/coordinates_4.txt ./01-bed_files/geneNames_4.txt ./01-bed_files/empty_score_4.txt ./01-bed_files/strand_4.txt > ./01-bed_files/nonSigGenes_141_rep1.bed
# Replicate 2
paste ./01-bed_files/chrom_5.txt ./01-bed_files/coordinates_5.txt ./01-bed_files/geneNames_5.txt ./01-bed_files/empty_score_5.txt ./01-bed_files/strand_5.txt > ./01-bed_files/nonSigGenes_27_rep2.bed
paste ./01-bed_files/chrom_6.txt ./01-bed_files/coordinates_6.txt ./01-bed_files/geneNames_6.txt ./01-bed_files/empty_score_6.txt ./01-bed_files/strand_6.txt > ./01-bed_files/nonSigGenes_8_rep2.bed
paste ./01-bed_files/chrom_7.txt ./01-bed_files/coordinates_7.txt ./01-bed_files/geneNames_7.txt ./01-bed_files/empty_score_7.txt ./01-bed_files/strand_7.txt > ./01-bed_files/nonSigGenes_172_rep2.bed
paste ./01-bed_files/chrom_8.txt ./01-bed_files/coordinates_8.txt ./01-bed_files/geneNames_8.txt ./01-bed_files/empty_score_8.txt ./01-bed_files/strand_8.txt > ./01-bed_files/nonSigGenes_141_rep2.bed
# Replicate 3
paste ./01-bed_files/chrom_9.txt ./01-bed_files/coordinates_9.txt ./01-bed_files/geneNames_9.txt ./01-bed_files/empty_score_9.txt ./01-bed_files/strand_9.txt > ./01-bed_files/nonSigGenes_27_rep3.bed
paste ./01-bed_files/chrom_10.txt ./01-bed_files/coordinates_10.txt ./01-bed_files/geneNames_10.txt ./01-bed_files/empty_score_10.txt ./01-bed_files/strand_10.txt > ./01-bed_files/nonSigGenes_8_rep3.bed
paste ./01-bed_files/chrom_11.txt ./01-bed_files/coordinates_11.txt ./01-bed_files/geneNames_11.txt ./01-bed_files/empty_score_11.txt ./01-bed_files/strand_11.txt > ./01-bed_files/nonSigGenes_172_rep3.bed
paste ./01-bed_files/chrom_12.txt ./01-bed_files/coordinates_12.txt ./01-bed_files/geneNames_12.txt ./01-bed_files/empty_score_12.txt ./01-bed_files/strand_12.txt > ./01-bed_files/nonSigGenes_141_rep3.bed

#
# FIFTH: DELETE UNNECESSARY INTERMEDIATE FILES
#####################################################################
rm ./01-bed_files/chrom_1.txt ./01-bed_files/chrom_2.txt ./01-bed_files/chrom_3.txt ./01-bed_files/chrom_4.txt 
rm ./01-bed_files/chrom_5.txt  ./01-bed_files/chrom_6.txt ./01-bed_files/chrom_7.txt ./01-bed_files/chrom_8.txt 
rm ./01-bed_files/chrom_9.txt ./01-bed_files/chrom_10.txt ./01-bed_files/chrom_11.txt ./01-bed_files/chrom_12.txt
rm ./01-bed_files/coordinates_1.txt ./01-bed_files/coordinates_2.txt ./01-bed_files/coordinates_3.txt ./01-bed_files/coordinates_4.txt 
rm ./01-bed_files/coordinates_5.txt ./01-bed_files/coordinates_6.txt ./01-bed_files/coordinates_7.txt ./01-bed_files/coordinates_8.txt
rm ./01-bed_files/coordinates_9.txt ./01-bed_files/coordinates_10.txt ./01-bed_files/coordinates_11.txt ./01-bed_files/coordinates_12.txt
rm ./01-bed_files/empty_score_1.txt ./01-bed_files/empty_score_2.txt ./01-bed_files/empty_score_3.txt ./01-bed_files/empty_score_4.txt 
rm ./01-bed_files/empty_score_5.txt ./01-bed_files/empty_score_6.txt ./01-bed_files/empty_score_7.txt ./01-bed_files/empty_score_8.txt
rm ./01-bed_files/empty_score_9.txt ./01-bed_files/empty_score_10.txt ./01-bed_files/empty_score_11.txt ./01-bed_files/empty_score_12.txt
rm ./01-bed_files/strand_1.txt ./01-bed_files/strand_2.txt ./01-bed_files/strand_3.txt ./01-bed_files/strand_4.txt
rm ./01-bed_files/strand_5.txt ./01-bed_files/strand_6.txt ./01-bed_files/strand_7.txt ./01-bed_files/strand_8.txt
rm ./01-bed_files/strand_9.txt ./01-bed_files/strand_10.txt ./01-bed_files/strand_11.txt ./01-bed_files/strand_12.txt
rm ./01-bed_files/IDinfo_1.txt ./01-bed_files/IDinfo_2.txt ./01-bed_files/IDinfo_3.txt ./01-bed_files/IDinfo_4.txt 
rm ./01-bed_files/IDinfo_5.txt ./01-bed_files/IDinfo_6.txt ./01-bed_files/IDinfo_7.txt ./01-bed_files/IDinfo_8.txt
rm ./01-bed_files/IDinfo_9.txt ./01-bed_files/IDinfo_10.txt ./01-bed_files/IDinfo_11.txt ./01-bed_files/IDinfo_12.txt
rm ./01-bed_files/geneNames_1.txt ./01-bed_files/geneNames_2.txt ./01-bed_files/geneNames_3.txt ./01-bed_files/geneNames_4.txt 
rm ./01-bed_files/geneNames_5.txt ./01-bed_files/geneNames_6.txt ./01-bed_files/geneNames_7.txt ./01-bed_files/geneNames_8.txt
rm ./01-bed_files/geneNames_9.txt ./01-bed_files/geneNames_10.txt ./01-bed_files/geneNames_11.txt ./01-bed_files/geneNames_12.txt
rm ./01-bed_files/nonSig_27DEG_rep1_results.txt ./01-bed_files/nonSig_27DEG_rep2_results.txt ./01-bed_files/nonSig_27DEG_rep3_results.txt
rm ./01-bed_files/nonSig_8DEG_rep1_results.txt ./01-bed_files/nonSig_8DEG_rep2_results.txt ./01-bed_files/nonSig_8DEG_rep3_results.txt
rm ./01-bed_files/nonSig_172DEG_rep1_results.txt ./01-bed_files/nonSig_172DEG_rep2_results.txt ./01-bed_files/nonSig_172DEG_rep3_results.txt
rm ./01-bed_files/nonSig_141DEG_rep1_results.txt ./01-bed_files/nonSig_141DEG_rep2_results.txt ./01-bed_files/nonSig_141DEG_rep3_results.txt    
#
# END
