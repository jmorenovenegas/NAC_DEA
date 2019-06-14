## DEA Core script 3. Differential Expression Analysis
# Input: matrix with normalised counts. ncounts
# Input: dataframe with samples to analyse 
# Output: table with differential expression data according to feature.

# Based on `https://github.com/hbc/clark-nanostring`, `https://rawgit.com/hbc/clark-nanostring/master/tumor-set/tumor-set.html#t-test`, `http://bioconductor.org/packages/devel/bioc/vignettes/NanoStringQCPro/inst/doc/vignetteNanoStringQCPro.pdf`

##### Libraries ######
loadpkg("NanoStringQCPro")
loadpkg("dplyr")
loadpkg("matrixStats")
loadpkg("tidyr")
loadpkg("ggplot2")
loadpkg("cowplot")
loadpkg("DESeq2")

getDE <-function(ncounts, metadata, design){
	dds = DESeqDataSetFromMatrix(countData = ncounts, colData = metadata, design = design)
  sizeFactors(dds) = rep(1, ncol(ncounts))
  dds = DESeq(dds)
  res = results(dds) %>%
        data.frame() %>%
        tibble::rownames_to_column() %>%
        arrange(pvalue) %>%
        mutate(rowname=gsub("NM_", "NM-", rowname)) %>%
        tidyr::separate(rowname, c("type", "gene.name", "id"), sep="_") %>%
        mutate(id=gsub("NM-", "NM_", id))
  return(res)
}

## Gets the same DE data, less processing, for GSEA in coreDEA4.R
getDE.raw <-function(ncounts, metadata, design){ 
	dds <-  DESeqDataSetFromMatrix(countData = ncounts, colData = metadata, design = design)
  sizeFactors(dds) <-  rep(1, ncol(ncounts))
  dds <-  DESeq(dds)
  res <-  results(dds) %>%
          data.frame() %>%
          tibble::rownames_to_column() %>% 
          arrange(pvalue)
  return(res)
}



