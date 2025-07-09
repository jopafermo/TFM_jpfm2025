#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: This script is used to find known motifs for EIN3 transcription factor 
# in the promoter sequences of the known unique upregulated genes from the heatmap analysis
# Bioinformatic tools: FUZZNUC (from EMBOSS: https://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html)
# File requirement: this script will uses the output from two previous scripts:
# - TFBS-analysis_PROMOTER_isolation.sh
# - TFBS-analysis_BED_Genome_&_RandmSet.sh
# Additionally a motif pattern is required to run the search:
# - EBSnew = [AG]T[ATG]CA[AT]TG[ATC]A[TC]
#///////////////////////////////////////////////////////////////////////
#
## 1@. SEARCHING ON PROMOTERS FROM THE DEG SETS
###########################################################################
# I'm using -rformat2 excel to get the output in a format that can be easily read in Excel and manipulated to extract information.
# removing the -rformat2 option, the output will be, by default, a seqtable, which is difficult to manipulate, but it has at the end a total count of sequences with #hits.
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_G0hg1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_G0hg2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_G48hg3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_G48hg4_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg4_filtered_EBSn.fuzznuc
#
## 2@. SEARCHING ON PROMOTERS FROM THE NOSIG-DEG SETS
###########################################################################
# Replicate 1
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_27_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_8_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_172_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_141_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep1_filtered_EBSn.fuzznuc
# Replicate 2
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_27_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_8_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_172_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_141_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep2_filtered_EBSn.fuzznuc
# Replicate 3
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_27_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_8_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_172_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_noSigDEG_141_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep3_filtered_EBSn.fuzznuc
#
## 3@. SEARCHING ON PROMOTERS FROM THE RANDOM SETS
###########################################################################
# Replicate 1
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_27random_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_8random_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_172random_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep1_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_141random_rep1_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep1_filtered_EBSn.fuzznuc
# Replicate 2
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_27random_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_8random_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_172random_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep2_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_141random_rep2_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep2_filtered_EBSn.fuzznuc
# Replicate 3
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_27random_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_8random_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_172random_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep3_filtered_EBSn.fuzznuc
fuzznuc -snucleotide_sequence ./02-promoters_files/2Kprom_141random_rep3_filtered.fa.out -pattern [AG]T[ATG]CA[AT]TG[ATC]A[TC] -complement -rformat2 excel -outfile ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep3_filtered_EBSn.fuzznuc
#
# END
