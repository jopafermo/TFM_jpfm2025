# TFM_jpfm2025_DGE-RNAseq_and_Promoter_Motif_Search
Scripts adapted, optimized and/or created from scratch used to complete the TFM

Folder 1: 01-dge-rnaseq_analysis 
  This folder contains the R-Notebook with the entire pipeline for DGE-RNAseq analysis performed in RStudio.

Folder 2: 02-genomic_pipelines_for_promoter_isolation
  This folder include the different scripts (bash) used to isolate the promoter sequence of candidate genes
  It also contains a README file (README_EXTRACCION_PROMOTORES.rtf) with additional information about the samples and the purpose of the scripts.

Folder 3: 03-targeted_EBS_motif_search_fuzznuc
  This folder harbors the different scripts (bash) used to scan the 2EBS(-1) motif in the isolated promoter regions using fuzznuc (EMBOSS).
  Similarly to the previous foulder, it also contains a README file (README_RESULTADOS_FUZZNUC.rtf) with additional information about the samples and the purpose of the scripts.
  
Folder 4: 04-exploratory _motif_search_meme
  Finally, this folder include the different scripts (bash) used to explore different motifs in the isolated promoter regions by an differential enrichment approach using meme (MEME SUITE).
  Like the other two precious folders, it also contains a README file (README_RESULTADOS_MEME.rtf) with additional information about the samples and the purpose of the scripts.

Each one of these scripts has a full description of each step within.

Josefina Patricia Fernández Moreno
09MBIF - TRABAJO FIN DE MASTER
VIU, 2024-2025
