## This is a differential expression analysis of Basal-like breast cancer tumour samples to compare gene expression patterns before and after de NAC (NeoAdyuvant Chemotherapy).
## The protocol is base on the execution of three R scripts that perform the following steps:
##  1. Load RCC files. Build counts matrix
##  2. Data normalisation
##  3. DE analysis


#######################################################
# SETTING UP R SESSION
#######################################################

# set seed to be reproduciblish
set.seed(12345)

# set working directory
workingDir <- "."
setwd(workingDir)

options(stringsAsFactors = FALSE)

# Load required libraries
source("Code/coreDEA_1.R")
source("Code/coreDEA_2.R")
source("Code/coreDEA_3.R")
source("Code/coreDEA_4.R")
loadpkg("rgl")
loadpkg("calibrate") # To label the volcano plot


#######################################################
# LOAD DATA
#######################################################

# Get a matrix with the raw counts from the set of samples we want to analyse. 
# In this case we select Basal-like NAC samples and their correspondent RCC files to build an RCCSet and get the counts matrix

DataDir <- paste("Data", sep = '/')
rccDir <- file.path(DataDir, "RCCData")

## Process file names: <sample>_<key1>_<key2>_<subtype>.RCC
## We have 204 files. There are 6 files without <subtype> : <sample>_<key1>_<key2>.RCC
## Key 1 contains pre/post info. 0 = PRE ; 1 = POST
metadata <-  data.frame(fname = list.files(path = rccDir, pattern="*.RCC")) %>%
  tidyr::separate(fname,c("Sample.name","key1","key2",'other'),sep="_",remove = F) %>%
  tidyr::separate(other, c('subtype', 'ext'), sep = '\\.', remove = T) %>%
  dplyr::filter((key1 == 0 | key1 == 1) & subtype == 'BL') %>% #Ignore Key2 (pCR info)
  dplyr::select(fname, Sample.name, key1, subtype) 

eset <- getRCCSet(metadata)

#######################################################
# QC (Quality Control)
#######################################################

# Even though the samples we analyse here have passed the QC checks via nSolver. We monitor some QC aspects and normalise the data according to several criteria.

##### Imaging QC
# Imaging QC refers to the percentage of fields of view (FOVs) successfully counted by a Digital Analyzer scan of a lane. When a substantial percentage of FOVs are not successfully counted, there may be issues with the resulting data (see nSolver User Manual).
pdf(file = 'Results/NAC-BL/QC_plotFOV.pdf')
plotFOV.key1(eset = eset, metadata = metadata, fov_threshold = 80)
dev.off()

# As there is a sample under the % FOV counted threshold we remove it
metadata <- pData(eset) %>% 
  tibble::rownames_to_column(var = "Sample.name") %>%
  inner_join(metadata, by="Sample.name") %>%
  dplyr::filter(FovCounted/FovCount * 100 > 80) %>%
  dplyr::select(fname, Sample.name, key1)


##### Binding density QC
#Due to the nature of the nCounter technology, analysis of some samples may produce too many or too few probes to be accurately counted by the Digital Analyzer. When too many probes are present, the Digital Analyzer is not able to distinguish each and every probe present in the lane. When too few fluorescent species are present, the Digital Analyzer may have difficulty focusing on the lane surface. Therefore, a measurement of mean binding density (spots per square micron) is provided with each lane scanned. The linear range of counting extends from 0.05 to 2.25 spots per square micron for assays run on an nCounter MAX or FLEX system. The range is 0.05 to 1.18 spots per square micron for assays run on the nCounter SPRINT system.  (see nSolver User Manual).
pdf(file = 'Results/NAC-BL/QC_plotBD.pdf')
plotBD.key1(eset = eset, metadata = metadata)
dev.off()
# All samples are in the binding density range


##### Positive Controls
# Check the expression of the positive controls. The positive genes follow the expeted pattern of expresssion.
pdf(file = 'Results/NAC-BL/QC_posContrl_boxplot.pdf')
boxplot_expr(eset,is_positive)
dev.off()

# Check the expression of the negative controls
pdf(file = 'Results/NAC-BL/QC_negContrl_boxplot.pdf')
boxplot_expr(eset,is_negative)
dev.off()

##### Noise Threshold  
# We establish a noise threshold. This threshold is based on the mean and standard deviation of counts of the negative control genes and represents the background noise. We define it as the mean expression of the negative genes counts + 2 times the standard deviation.
lodcounts = extract_pred(eset, is_negative)
lod = mean(lodcounts$count) + 2 * sd(lodcounts$count)


##### Housekeeping genes  
# Expression of each housekeeping genes in all samples. The line in red represents the noise threshold.
pdf(file = 'Results/NAC-BL/hK_exprs_boxplot.pdf')
boxplot_expr(eset,is_housekeeping) + geom_hline(yintercept = (lod),colour="red")
dev.off()
# We can observe that all the housekeeping genes have enough expression (above the noise threshold).


##### Expression of all the housekeeping genes in each sample.  
# We plot the mean expression of all the housekeeping genes in each sample.
pdf(file = 'Results/NAC-BL/hK_mean_exprs_boxplot.pdf')
hK_mean_exprs_plot(eset, metadata, lod)
dev.off()
# We can obverve that some sample have overall low expression for the housekeeping genes, but all are above the noise threshold.


##### General expression
# Expression of all genes (Endogenous + Housekeeping).
pdf(file = 'Results/NAC-BL/general_exprs_plot.pdf')
general_exprs_plot(eset, metadata, noise_threshold = lod)
dev.off()


#######################################################
# NORMALIZATION
#######################################################

##### Housekeeping normalization
# We select thouse housekeeping genes that have expression values greater than the noise threshold and a mean value of expression of at least 200 counts.

hk = counts[grepl("Housekeeping", rownames(counts)),]
abovenoise = rowSums(hk > (lod)) >= (ncol(hk))
hk_abovenoise = hk[abovenoise,]
aboveMean = (apply(hk_abovenoise,1,mean))>= 200
hk_sel= hk_abovenoise[aboveMean,]
hk_norm=rownames(hk_sel)

hk_pos = function(counts, hk_norm) {
  hk = counts[hk_norm,]
  geoMeans = apply(hk, 2, function(col) exp(mean(log(col[col != 0]))))
  return(geoMeans)}
hkFactor = function(counts, hk_norm) {
  geoMeans = hk_pos(counts, hk_norm)
  nf = mean(geoMeans) / geoMeans
  return(nf)}

metadata$hk_nf = hkFactor(counts, hk_norm)
ncounts = counts %*% diag(metadata$hk_nf)
colnames(ncounts) = colnames(counts)
postnorm = ncounts %>%
  data.frame() %>%
  tidyr::gather("sample", "count")
postnorm$sample<-gsub("X","",postnorm$sample)

pdf(file = 'Results/NAC-BL/postnorm_exprs_plot.pdf')
ggplot(postnorm, aes(sample, count)) + geom_boxplot(colour = "black", fill = "#56B4E9",outlier.size = 0.5) +
  scale_y_continuous(trans = "log2") + 
  xlab("") + ggtitle("") + theme(axis.text.x = element_text(angle = 90, 
                                                            hjust = 1,vjust = 0.5,size = 5))+
  geom_smooth(se=T, aes(group=1))+
  geom_hline(yintercept = lod,colour="red")
dev.off()

##### Drop genes and produce the normalised counts matrix
# We’ll drop the genes that are below the LOD in over 80% of the samples:
allNames <- rownames(ncounts)
ncounts <-  ncounts[(rowSums(ncounts < lod) < round((0.2 * ncol(ncounts)),0)),]
filteredNames <- rownames(ncounts)
filterOutGenes <- allNames[allNames%!in%filteredNames]

# We save the metadata dataframe and the normalised counts matrix for further analysis.
save(metadata, file = "Data/RData/BL-NAC_metadata.RData")
save(ncounts, file = "Data/RData/BL-NAC_HKnormalised_counts.RData")

#######################################################
# DEA
####################################################### 

# DEA by the DESeq2 library. The ncounts matrix is converted to integers for the analysis. The original counts are preserved in counts(dds). 
# See [this thread](https://www.biostars.org/p/101727/) for log2FC.

load(file = "Data/RData/BL-NAC_metadata.RData")
load(file = "Data/RData/BL-NAC_HKnormalised_counts.RData")

design = ~ key1
dea.o <- getDE(ncounts = round(ncounts), metadata = metadata, design = design)
write_csv(dea.o, "Results/NAC-DEA.csv")
dea.p <- dea.o %>% dplyr::filter(padj <0.05) %>% arrange(-log2FoldChange)
View(dea.p)

# We can show the results of the differential expression analysis as a volcano plot. 
# A volcano plot typically plots some measure of effect on the x-axis (typically the fold change) and the statistical significance on the y-axis (typically the -log10 of the p-value). 
# Genes that are highly dysregulated are farther to the left and right sides, while highly significant changes appear higher on the plot.

tab = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj))
tab2 = data.frame(logFC = dea.o$log2FoldChange, negLogPval = -log10(dea.o$padj), Gene=dea.o$gene.name)
lfc = 1
pval = 0.05

pdf(file = 'Results/NAC-BL/volcano_plot.pdf')
par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))

points(tab[(abs(tab$logFC) > lfc), ], pch = 16, cex = 0.8, col = "orange") 
points(tab[(tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "green") 
points(tab[(abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "red") 
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

with(subset(tab2, negLogPval > -log10(pval) & abs(logFC)>1), textxy(logFC, negLogPval, labs=Gene, cex=.4))
dev.off()

# Save significant genes into a file
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
write.csv(tab2[which(signGenes),], file = "Results/NAC-BL/signGenes.csv")


#######################################################
# PATHWAY ANALYSIS
#######################################################

# Here we use GSEA to perform pathway analysis on the list of genes. 
# Gene Set Enrichment Analysis focuses on the changes of expression in groups of genes, and by doing so, this method resolves the problem of the undetectable, small changes in the expression of single genes [see wikipedia.](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis).

load(file = "Data/RData/BL-NAC_metadata.RData")
load(file = "Data/RData/BL-NAC_HKnormalised_counts.RData")

design = ~ key1
res <- getDE.raw(ncounts = round(ncounts), metadata = metadata, design = design) # the core script for the pathway analysis expects a dataframe with the DEA named res
converted = unlist(lapply(strsplit(res$rowname, "_", fixed = TRUE), "[", 4))
converted = unlist(lapply(strsplit(converted, ".", fixed = TRUE), "[", 1))
converted = paste0("NM_", converted)
res$refseq_mrna = converted

lost.genes = length(res$rowname) - length((res %>% data.frame() %>% left_join(entrez, by = "refseq_mrna") %>% 
                                             filter(!is.na(entrezgene)) %>% filter(!is.na(log2FoldChange)) %>% filter(!is.na(lfcSE)))$rowname)

gsea_rs = gsea_cp(res, "NACstatus")

gsea_summary = gsea_rs$summary %>% arrange(pvalue)
gsea_summary = convert_gsea_results(gsea_summary)
write_csv(x = gsea_summary, path = "Results/NAC-BL/BL-NACstatus_gsea.csv")

# We need to first pull out the Refseq IDs from the Nanostring IDs and convert those to Entrez IDs. Not all of the genes have Entrez IDs. 
# Genes are sorted by signal to noise ratio, which are defined as the log2Fold change (the signal) divided by the log2Fold change standard error (the noise).
# I sorted the ontology/KEGG pathway terms by p-value and output a table of the results. 
# The core_symbols are the genes that most strongly contributed to the pathway score, so the ones with the highest or lowest enrichment score depending on the direction of enrichment.

# We perform a classic enrichment analysis in GO, Reactome and KEGG for all the genes differentially epressed by NAC.
enrich_rs = enrich_cp(res, "NAC", type="all")
enrich_summary = enrich_rs$summary %>% arrange(p.adjust)
enrich_summary = convert_enriched_ids(enrich_summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrich_summary,path = "Results/NAC-BL/BL-NACstatus_enrichment.csv" )

# We perform also a classic enrichment analysis for GO, Reactome and KEGG using the genes substantially over-expressed and under-expressed. We will use a two-fold fold-change cutoff
enrichover_rs = enrich_cp(res, "NAC", type="over")
enrichover_summary = enrichover_rs$summary %>% arrange(p.adjust)
enrichover_summary = convert_enriched_ids(enrichover_summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrichover_summary,path = "Results/NAC-BL/BL-NACstatus_over_enrichment.csv" )

enrichunder_rs = enrich_cp(res, "NAC", type="under")
enrichunder_summary = enrichunder_rs$summary %>% arrange(p.adjust)
enrichunder_summary = convert_enriched_ids(enrichunder_summary,entrezsymbol = entrezsymbol) %>% arrange(p.adjust)
write_csv(x = enrichunder_summary,path = "Results/NAC-BL/BL-NACstatus_under_enrichment.csv" )

# Finally we can see the GO BP, Reactome and KEGG pathways enriched.
save(enrich_rs, enrichunder_rs, enrichover_rs, file = 'Data/NAC-BL_enrichments.RData')
load(file = 'Data/NAC-BL_enrichments.RData')

# GO BP Enrichment in all the significantly expressed genes:
pdf(file = 'Results/NAC-BL/GO_BP_all_dotplot.pdf')
dotplot(enrich_rs$bp, x="count", showCategory=10, color="qvalue")
dev.off()

# Reactome enrichment in all significantly expressed genes:
pdf(file = 'Results/NAC-BL/Reactome_all_dotplot.pdf')
dotplot(enrich_rs$re, x="count", showCategory=10, color="qvalue")
dev.off()

# KEGG enrichment in all significantly expressed genes:
pdf(file = 'Results/NAC-BL/KEGG_all_dotplot.pdf')
dotplot(enrich_rs$kg, x="count", showCategory=10, color="qvalue")
dev.off()

# GO BP Enrichment in over-expressed genes:
pdf(file = 'Results/NAC-BL/GO_BP_over_dotplot.pdf')
dotplot(enrichover_rs$bp, x="count", showCategory=10, color="qvalue")
dev.off()


# Reactome enrichment in over-expressed genes:
pdf(file = 'Results/NAC-BL/Reactome_over_dotplot.pdf')
dotplot(enrichover_rs$re, x="count", showCategory=10, color="qvalue")
dev.off()

# KEGG enrichment in over-expressed genes:
pdf(file = 'Results/NAC-BL/KEGG_over_dotplot.pdf')
dotplot(enrichover_rs$kg, x="count", showCategory=10, color="qvalue")
dev.off()


# GO BP Enrichment in under-expressed genes:
pdf(file = 'Results/NAC-BL/GO_BP_under_dotplot.pdf')
dotplot(enrichunder_rs$bp, x="count", showCategory=10, color="qvalue")
dev.off()


# Reactome enrichment in under-expressed genes:
pdf(file = 'Results/NAC-BL/Reactome_under_dotplot.pdf')
dotplot(enrichunder_rs$re, x="count", showCategory=10, color="qvalue")
dev.off()

# KEGG enrichment in under-expressed genes:
pdf(file = 'Results/NAC-BL/KEGG_under_dotplot.pdf')
dotplot(enrichunder_rs$kg, x="count", showCategory=10, color="qvalue")
dev.off()