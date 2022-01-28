# BIN2_RAPTOR1B_MULTIOMICS
Code repository for Montes et al., 2022

Brief description of files:

`BIN2_MAKS_Phospho_TMT_normalization_and_DE_edgeR.Rmd`  
R markdown script for BIN2 MAKS phosphoproteomics data analysis which performs:  
- Sample loading (SL), and trimmed mean of M (TMM) normalization
- Differential phosphosite intensities assessment using edgeR  
The script will output a html file including many diagnostic plots  

`BIN2_RAPTOR_mutants_Phospho_TMT_IRS_normalization_and_DE_edgeR.Rmd`  
R markdown script for bin2D, bin2T, and raptor1b mutants phosphoproteomics data analysis which performs:  
- Sample loading (SL), trimmed mean of M (TMM), and internal reference scaling (IRS) normalization
- Differential phosphosite intensities assessment using edgeR  
The script will output a html file including many diagnostic plots  

`Fig.5a_heatmap.R`  
R code used to plot Fig. 5 heatmap. The script will output a PDF file  

`NMC_phenotyping_MLM.R`  
Script used for BL phenotype and MDC staining under sucrose starvation statistical testing by generalized linear model regression

`QuantSeq_dataprocess.slurm`  
Script for SLURM manager that was used to perform the QuantSeq data processing:  
- Trimming of raw reads using BBDUK
- Reads quality check using fastQC  
- Sequence alignment to TAIR10 annotation using STARaligner
- Indexing the alignment files using Samtools
- Extract transcript counts using HTSeq-count

`Transcriptomic_analysis_PoissonSeq.Rmd`  
R markdown script for bin2D, bin2T, and raptor1b mutants transcript differential expression analysis:
- Use TMM normalization
- Use PoissonSeq for differential expression  
The script will output a html file including many diagnostic plots

