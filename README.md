README.md

In this repository differential expression analysis(DEA) and pathway analysis are done over gene expression data from breast cancer patients samples. 
The expression is cuantified using NanoString's nCounter technology.

In Data/RCCData directory there are the .RCC files with gene expression information corresponding to each sample.

In Code/ directory there are R scripts with functions used in main R scripts.

To perform the analysis over all cancer subtypes run NAC_status_DEA.R. The RData files will be saved in Data/RData directory and the results will be saved in Results/NAC/ directory.

To perform the analysis only over Basal-like cancer subtype run NAC-BL_status_DEA.R. The RData files will be saved in Data/RData directory and the results will be saved in Results/NAC-BL/ directory.

The process is the same in both analysis. First, the data pass a quality control and it is normalized. After that DEA is done and finally pathway analysis is done.

After the execution you could find:
- a csv file with DEA analysis results.
- a volcano plot
- a csv file with significant genes
- 2 RData files with metadata and Housekeeping normalised counts
- a csv file with gene set enrichment analysis
- 3 csv files with summary enrichments (all, under-expressed and over-expressed)
- 3 GO BP enrichment dotplots (all, under-exressed, over-expressed)
- 3 KEGG enrichment dotplots (all, under-exressed, over-expressed)
- 3 Reactome enrichment dotplots (all, under-exressed, over-expressed)

Requirements:
- R version 3.5.1
R libraries (installations calls inside the code):
- NanoStringQCPro
- dplyr
- readr
- matrixStats
- tidyr
- ggplot2
- cowplot
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- biomaRt
- pathview
- ReactomePA
- DOSE
- rgl
- calibrate
