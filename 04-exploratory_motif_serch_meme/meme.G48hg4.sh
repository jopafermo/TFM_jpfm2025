#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025, June
# Goal: This script is used to explore binding site motifs with MEME in promoter regions
# Bioinformatic tools: MEME (https://meme-suite.org/meme/doc/meme.html?man_type=web), from MEME SUITE (https://meme-suite.org/meme/index.html)
# File required: the fasta files (.fa.out) with the promoter sequences obtained from the genomics pipelines for G48hg4 and the six control files for 141 genes.
#///////////////////////////////////////////////////////////////////////
#
## 1@. DIFFERENTIAL ENRICHMENT ANALYSIS
##########################################################################################
# With the mode "se", meme will perform a search of motifs in both DEG set G48hg4 and each one of the control set replicates, 
# then will contrast the results and it will provide those motifs enriched only in G48hg4.
#
# DEG-G48hg4 vs NonSigDEG141 (rep1-rep3) controls
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_noSigDEG_141_rep1_filtered.fa.out -o noSigDEG141_rep1 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_noSigDEG_141_rep2_filtered.fa.out -o noSigDEG141_rep2 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_noSigDEG_141_rep3_filtered.fa.out -o noSigDEG141_rep3 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
# DEG-G48hg4 vs random141 (rep1-rep3) controls
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_141random_rep1_filtered.fa.out -o random141_rep1 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_141random_rep2_filtered.fa.out -o random141_rep2 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
meme ./2Kprom_G48hg4_filtered.fa.out -neg ./2Kprom_141random_rep3_filtered.fa.out -o random141_rep3 -dna -nostatus -time 1440 -mod anr -nmotifs 30 -minw 6 -maxw 12 -objfun se -minsites 1 -revcomp -markov_order 0
#
# Description of the elements of the code:
# -neg = contol sequence
# -o = output fule
# -dna = specification of the alphabet for DNA sequences
# -nostatus = print no status messages to terminal
# -time = maximum time invested in completing the code (by default is 1440 seconds)
# -mod anr = assumption of any number of repetitions when finding hits
# -nmotifs =  maximum number of motifs to find
# -minw = minimum size of the motif (in base pairs)
# -maxw = maximum size of the motif (in base pairs)
# -objfun se = selective differential enrichment
# -minsites 1 =  find motifs with at least 1 hit
# -markov_order 0 = is the default value used by meme as maximum order of the Markov model to read from background files
#
# END