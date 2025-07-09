#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: The goal of this script is to obtain the coordinates of both gene and promoters for all the genes 
# within the Medicago truncatula genome and to create randomize sets of promoters from the entire genome.
# Bioinformatic tools: bedtools (https://bedtools.readthedocs.io/en/latest/content/overview.html)
# File requirement: this script will uses the output from two previous scripts:
# - TFBS-analysis_BED_formatting.sh
# - TFBS-analysis_BED_Genome_&_RandmSet.sh
#///////////////////////////////////////////////////////////////////////
#
## 1@. OBTAINING THE COORDINATES OF THE PROMOTERS FOR THE DEG SET
##########################################################
# FIRST, OBTAINING THE COORDINATES OF THE PROMOTERS
# I use bedtools flank to obtain the coordinates for the 2kb of the promoters with the option -s to choose the strand accordingly
bedtools flank -i ./01-bed_files/A17GeneUpG0hg1.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_G0hg1.bed
bedtools flank -i ./01-bed_files/A17GeneUpG0hg2.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_G0hg2.bed
bedtools flank -i ./01-bed_files/A17GeneUpG48hg3.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_G48hg3.bed
bedtools flank -i ./01-bed_files/A17GeneUpG48hg4.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_G48hg4.bed
#
# SECOND, CUTTING OVERLAPPING PROMOTER SEQUENCES WITH SURROUNDING GENE SEQUENCES
# I have to extract the coordinates. For that I use bedtools subtract with the -s option to keep the strand information
bedtools subtract -a ./02-promoters_files/2Kprom_G0hg1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_G0hg1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_G0hg2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_G0hg2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_G48hg3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_G48hg3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_G48hg4.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_G48hg4_nooverlap.bed 
# It is possible that some of the promoters are completely overlapped by the gene sequences, so I have to filter them by a minimum length
#
# THIRD, FILTERING PROMOTERS BY LENGTH
# I will use awk to filter the promoters by length, keeping only those with a length equal or greater than 600nt
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_G0hg1_nooverlap.bed > ./02-promoters_files/2Kprom_G0hg1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_G0hg2_nooverlap.bed > ./02-promoters_files/2Kprom_G0hg2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_G48hg3_nooverlap.bed > ./02-promoters_files/2Kprom_G48hg3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_G48hg4_nooverlap.bed > ./02-promoters_files/2Kprom_G48hg4_filtered.bed
#
# FOURTH, OBTAINING THE NUCLEOTIDE SEQUENCE FOR THE FILTERED PROMOTER LIST
# I will obtain the nucleotide sequence of the filtered promoters using bedtools getfasta
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_G0hg1_filtered.bed -fo ./02-promoters_files/2Kprom_G0hg1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_G0hg2_filtered.bed -fo ./02-promoters_files/2Kprom_G0hg2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_G48hg3_filtered.bed -fo ./02-promoters_files/2Kprom_G48hg3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_G48hg4_filtered.bed -fo ./02-promoters_files/2Kprom_G48hg4_filtered.fa.out -name
# These files will be used to perform the motif searching with FUZZNUC
#
## 2@. OBTAINING THE COORDINATES OF THE PROMOTERS FOR THE NOSIG-DEG SET
##########################################################
# FIRST, OBTAINING THE COORDINATES OF THE PROMOTERS
# Replicate 1
bedtools flank -i ./01-bed_files/nonSigGenes_27_rep1.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep1.bed
bedtools flank -i ./01-bed_files/nonSigGenes_8_rep1.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep1.bed
bedtools flank -i ./01-bed_files/nonSigGenes_172_rep1.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep1.bed
bedtools flank -i ./01-bed_files/nonSigGenes_141_rep1.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep1.bed
# Replicate 2
bedtools flank -i ./01-bed_files/nonSigGenes_27_rep2.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep2.bed
bedtools flank -i ./01-bed_files/nonSigGenes_8_rep2.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep2.bed
bedtools flank -i ./01-bed_files/nonSigGenes_172_rep2.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep2.bed
bedtools flank -i ./01-bed_files/nonSigGenes_141_rep2.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep2.bed
# Replicate 3
bedtools flank -i ./01-bed_files/nonSigGenes_27_rep3.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep3.bed
bedtools flank -i ./01-bed_files/nonSigGenes_8_rep3.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep3.bed
bedtools flank -i ./01-bed_files/nonSigGenes_172_rep3.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep3.bed
bedtools flank -i ./01-bed_files/nonSigGenes_141_rep3.bed -g ./02-promoters_files/MtGenome.chromsizes -l 2000 -r 0 -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep3.bed
#
# SECOND, CUTTING OVERLAPPING PROMOTER SEQUENCES WITH SURROUNDING GENE SEQUENCES
# Replicate 1
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_27_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_8_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_172_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep1_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_141_rep1.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep1_nooverlap.bed 
# Replicate 2
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_27_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_8_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_172_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep2_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_141_rep2.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep2_nooverlap.bed 
# Replicate 3
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_27_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_27_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_8_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_8_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_172_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_172_rep3_nooverlap.bed 
bedtools subtract -a ./02-promoters_files/2Kprom_noSigDEG_141_rep3.bed -b ./01-bed_files/MtAllGenes.bed -s > ./02-promoters_files/2Kprom_noSigDEG_141_rep3_nooverlap.bed 
#
# THIRD, FILTERING PROMOTERS BY LENGTH
# Replicate 1
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_27_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_27_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_8_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_8_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_172_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_172_rep1_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_141_rep1_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_141_rep1_filtered.bed
# Replicate 
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_27_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_27_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_8_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_8_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_172_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_172_rep2_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_141_rep2_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_141_rep2_filtered.bed
# Replicate 3
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_27_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_27_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_8_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_8_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_172_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_172_rep3_filtered.bed
awk 'BEGIN{OFS="\t"} {if($3-$2 >= 600) print $0}' ./02-promoters_files/2Kprom_noSigDEG_141_rep3_nooverlap.bed > ./02-promoters_files/2Kprom_noSigDEG_141_rep3_filtered.bed
#
# FOURTH, OBTAINING THE NUCLEOTIDE SEQUENCE FOR THE FILTERED PROMOTER LIST
# Replicate 1
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_27_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_27_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_8_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_8_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_172_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_172_rep1_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_141_rep1_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_141_rep1_filtered.fa.out -name
# Replicate 2
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_27_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_27_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_8_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_8_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_172_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_172_rep2_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_141_rep2_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_141_rep2_filtered.fa.out -name
# Replicate 3
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_27_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_27_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_8_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_8_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_172_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_172_rep3_filtered.fa.out -name
bedtools getfasta -fi ./00-Selected_gene_list_files/MtrunA17r5.0-20161119-ANR.genome.fasta -bed ./02-promoters_files/2Kprom_noSigDEG_141_rep3_filtered.bed -fo ./02-promoters_files/2Kprom_noSigDEG_141_rep3_filtered.fa.out -name
#
# 3@. REMOVING UNNECESSARY INTERMEDIATE FILES
###########################################################
rm ./02-promoters_files/*_nooverlap.bed
rm ./02-promoters_files/*_filtered.bed

## END

