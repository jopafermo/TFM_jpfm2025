#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: This script is used to obtain the promoter sequences of the known unique upregulated genes from the heatmap analysis
# Bioinformatic tools: bedtools (https://bedtools.readthedocs.io/en/latest/content/overview.html)
# File requirement: this script will uses the output from the previous script
# - TFBS-analysis_DATA_prep.sh
#///////////////////////////////////////////////////////////////////////
#
## 1@. PREPARING THE GENOME.CHROMSIZES FILE (COORDINATES FOR EACH CHROMOSOME)
##################################################################################
# Selecting the headings of the chromosomes from the Medicago genome fasta file
grep ">" ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta > ./02-promoters_files/chrom_heading.txt
# Removing unnecessary formatting and transforming spaces into tabs so I can select columns
sed 's/>//g' < ./02-promoters_files/chrom_heading.txt | sed 's/ /\t/g' | sed 's/len=//g' > ./02-promoters_files/chrom_heading_columns.txt 
# Selecting the two first columns (chromosome name and length) and saving it as MtGenome.chromsizes
cut -f1,2 ./02-promoters_files/chrom_heading_columns.txt > ./02-promoters_files/MtGenome.chromsizes
# This file will be used as the genome size file for bedtools
#
## 2@. CREATING A BED FILE FOR GENE COORDINATES FROM THE WHOLE GENOME
##################################################################################
# In order to cut the overlapping sequences, I  have to create a bed file with the coordinates of the whole genome
# For that, I use gene.gff3 file, and I process it as I did for obtaining the coordinates of the knoun upregulated genes
cut -f1 ./01-bed_files/genes.gff3 > ./01-bed_files/chrom.txt
cut -f3,4 ./01-bed_files/genes.gff3  > ./01-bed_files/coordinates.txt
cut -f5 ./01-bed_files/genes.gff3 > ./01-bed_files/empty_score.txt 
cut -f6 ./01-bed_files/genes.gff3 > ./01-bed_files/strand.txt
cut -f7 ./01-bed_files/genes.gff3 > ./01-bed_files/IDinfo.txt
sed 's/;/\t/g' < ./01-bed_files/IDinfo.txt | cut -f1 | sed 's/:/\t/g' | cut -f2 > ./01-bed_files/geneNames.txt
paste ./01-bed_files/chrom.txt ./01-bed_files/coordinates.txt ./01-bed_files/geneNames.txt ./01-bed_files/empty_score.txt ./01-bed_files/strand.txt > ./01-bed_files/output.txt
tail -n +2 ./01-bed_files/output.txt > ./01-bed_files/MtAllGenes.bed # to remove the header line
#
## 3@. ISOLATING THE PROMOTER COORDINATES FROM THE MtAllGenes.bed FILE, CUTTING OVERLAPPING PROMOTER SEQUENCES & FILTERING BY LENGTH
##################################################################################
bedtools flank -i ./01-bed_files/MtAllGenes.bed  -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_ALL.bed
bedtools subtract -a ./02-promoters_files/2Kprom_ALL.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_ALL_nooverlap.bed 
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_ALL_nooverlap.bed > ./02-promoters_files/2Kprom_ALL_filtered.bed
#
## 4@. CREATING FILES WITH RANDOMLY SELECTED PROMOTER SEQUENCES WITH SIMILAR SIZES TO DEG GENES BY TRIPLICATE
###################################################################################
# Replicate 1
shuf -n 29 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_27random_rep1.bed
shuf -n 8 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_8random_rep1.bed
shuf -n 165 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_172random_rep1.bed
shuf -n 136 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_141random_rep1.bed
# Replicate 2
shuf -n 29 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_27random_rep2.bed
shuf -n 8 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_8random_rep2.bed
shuf -n 165 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_172random_rep2.bed
shuf -n 136 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_141random_rep2.bed
# Replicate 3
shuf -n 29 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_27random_rep3.bed
shuf -n 8 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_8random_rep3.bed
shuf -n 165 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_172random_rep3.bed
shuf -n 136 ./02-promoters_files/2Kprom_ALL_filtered.bed > ./02-promoters_files/2Kprom_141random_rep3.bed
#
## 5@. CUTTING OVERLAPPING PROMOTER SEQUENCES WITH GENE SEQUENCES FROM THE RANDOMLY SELECTED PROMOTERS
##################################################################################
# The promoter sequence that overlaps with an upstream gene sequence will be cut from the 2k promoter length
# Replicate 1
bedtools subtract -a ./02-promoters_files/2Kprom_27random_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_27random_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_8random_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_8random_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_172random_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_172random_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_141random_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_141random_rep1_nooverlap.bed 
# Replicate 2
bedtools subtract -a ./02-promoters_files/2Kprom_27random_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_27random_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_8random_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_8random_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_172random_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_172random_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_141random_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_141random_rep2_nooverlap.bed 
# Replicate 3
bedtools subtract -a ./02-promoters_files/2Kprom_27random_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_27random_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_8random_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_8random_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_172random_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_172random_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_141random_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_141random_rep3_nooverlap.bed 
#
## 6@. FILTERING THOSE PROMOTERS WITH SIZES SHORTER THAN 600 PAIR BASES
##################################################################################
# The non-overlapping promoters must be filtered by length, keeping only those with a length equal or greater than 600pb
# Replicate 1
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_27random_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_27random_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_8random_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_8random_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_172random_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_172random_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_141random_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_141random_rep1_filtered.bed
# Replicate 2
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_27random_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_27random_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_8random_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_8random_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_172random_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_172random_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_141random_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_141random_rep2_filtered.bed
# Replicate 3
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_27random_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_27random_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_8random_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_8random_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_172random_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_172random_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_141random_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_141random_rep3_filtered.bed
# 
## 7@. OBTAINING THE NUCLEOTIDE SEQUENCE FOR THE FILTERED RANDOM PROMOTERS
##########################################################
# I will obtain the nucleotide sequence of the filtered random promoters using bedtools getfasta
# Replicate 1
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_27random_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_27random_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_8random_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_8random_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_172random_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_172random_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_141random_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_141random_rep1_filtered.fa.out -name
# Replicate 2
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_27random_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_27random_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_8random_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_8random_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_172random_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_172random_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_141random_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_141random_rep2_filtered.fa.out -name
# Replicate 3
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_27random_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_27random_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_8random_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_8random_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_172random_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_172random_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_141random_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_141random_rep3_filtered.fa.out -name
# These files will be used to perform the motif searching with FUZZNUC as NEGATIVE CONTROL
#
# 10@. REMOVING UNNECESSARY INTERMEDIATE FILES
###########################################################
rm ./02-promoters_files/chrom_heading.txt
rm ./02-promoters_files/chrom_heading_columns.txt
rm ./01-bed_files/chrom.txt
rm ./01-bed_files/coordinates.txt
rm ./01-bed_files/empty_score.txt
rm ./01-bed_files/strand.txt
rm ./01-bed_files/IDinfo.txt
rm ./01-bed_files/geneNames.txt
rm ./01-bed_files/output.txt
rm ./02-promoters_files/*_nooverlap.bed
rm ./02-promoters_files/*_filtered.bed
#
## END