#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: xsss
# File requirement: outputs from TFBS-analysis_EBSnew-Motif_search.sh
#///////////////////////////////////////////////////////////////////////
#
## 1@. HEADINGS
########################################################################################################################
echo "Gene_Set" > ./04-EBS_hits_counting_files/col1.txt
echo "Init.tsv" > ./04-EBS_hits_counting_files/col2.txt 
echo "Filt.bed" > ./04-EBS_hits_counting_files/col3.txt 
echo "seq.EBSn" > ./04-EBS_hits_counting_files/col4.txt 
echo "EBSn.hit" > ./04-EBS_hits_counting_files/col5.txt 
paste ./04-EBS_hits_counting_files/col1.txt ./04-EBS_hits_counting_files/col2.txt ./04-EBS_hits_counting_files/col3.txt ./04-EBS_hits_counting_files/col4.txt ./04-EBS_hits_counting_files/col5.txt > ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/col*
#
## 2@. DEG DATASET
########################################################################################################################
echo "G0hg1" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/DEG.down.G0h.III.hmtree.genes_1.tsv | grep "downG" | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_G0hg1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "G0hg2" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/DEG.down.G0h.III.hmtree.genes_2.tsv | grep "downG" | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_G0hg2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G0hg2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "G48hg3" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/DEG.down.G48h.III.hmtree.genes_3.tsv | grep "downG" | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_G48hg3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm  ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "G48hg4" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/DEG.down.G48h.III.hmtree.genes_4.tsv | grep "downG" | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_G48hg4_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg4_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_G48hg4_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
## 2@. NOSIG-DEG DATASET
########################################################################################################################
# Replicate 1
echo "NoSigDEG_8_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_8DEG_rep1.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_8_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_27_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_27DEG_rep1.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_27_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_141_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_141DEG_rep1.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_141_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_172_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_172DEG_rep1.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_172_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
# Replicate 2
echo "NoSigDEG_8_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_8DEG_rep2.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_8_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_27_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_27DEG_rep2.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_27_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_141_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_141DEG_rep2.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_141_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_172_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_172DEG_rep2.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_172_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
# Replicate 3
echo "NoSigDEG_8_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_8DEG_rep3.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_8_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_8_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_27_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_27DEG_rep3.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_27_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_27_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_141_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_141DEG_rep3.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_141_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_141_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "NoSigDEG_172_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./00-Selected_gene_list_files/nonSig_172DEG_rep3.txt | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_noSigDEG_172_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_noSigDEG_172_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
## 3@. RANDOM DATASET
########################################################################################################################
# Replicate 1
echo "Random_8_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_8random_rep1.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_8random_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_27_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_27random_rep1.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_27random_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_141_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_141random_rep1.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_141random_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_172_rep1" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_172random_rep1.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_172random_rep1_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep1_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep1_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
# Replicate 2
echo "Random_8_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_8random_rep2.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_8random_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_27_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_27random_rep2.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_27random_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_141_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_141random_rep2.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_141random_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_172_rep2" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_172random_rep2.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_172random_rep2_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep2_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep2_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
# Replicate 3
echo "Random_8_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_8random_rep3.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_8random_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_8random_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_27_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_27random_rep3.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_27random_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_27random_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_141_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_141random_rep3.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_141random_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_141random_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
echo "Random_172_rep3" > ./04-EBS_hits_counting_files/N.txt
cat ./02-promoters_files/2Kprom_172random_rep3.bed | wc -l > ./04-EBS_hits_counting_files/A.txt
cat ./02-promoters_files/2Kprom_172random_rep3_filtered.fa.out | grep ">" | wc -l > ./04-EBS_hits_counting_files/B.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep3_filtered_EBSn.fuzznuc | grep "SeqName" | wc -l > ./04-EBS_hits_counting_files/C.txt
cat ./03-FUZZNUC_EBS_motive_files/2Kprom_172random_rep3_filtered_EBSn.fuzznuc | grep ":" | wc -l > ./04-EBS_hits_counting_files/D.txt
paste ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt >> ./04-EBS_hits_counting_files/EBShits_results.txt
rm ./04-EBS_hits_counting_files/N.txt ./04-EBS_hits_counting_files/A.txt ./04-EBS_hits_counting_files/B.txt ./04-EBS_hits_counting_files/C.txt ./04-EBS_hits_counting_files/D.txt
#
## END