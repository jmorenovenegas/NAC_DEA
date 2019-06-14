## DEA Core script 2. data normalisation
# Input: matrix with raw counts. counts
# Input: dataframe with samples to analyse 
# Output: This provides several functions to perform QC, plots and normalisation. The normalised counts matrix will be produced in the protocol script.

# Based on `https://github.com/hbc/clark-nanostring`, `https://rawgit.com/hbc/clark-nanostring/master/tumor-set/tumor-set.html#t-test`, `http://bioconductor.org/packages/devel/bioc/vignettes/NanoStringQCPro/inst/doc/vignetteNanoStringQCPro.pdf`

##### Libraries ######
source("0_loadlibraries.R")
loadpkg("NanoStringQCPro")
loadpkg("dplyr")
loadpkg("matrixStats")
loadpkg("tidyr")
loadpkg("ggplot2")
loadpkg("cowplot")
loadpkg("DESeq2")

#### Functions #####

'%!in%' <- function(x,y)!('%in%'(x,y))

pca_loadings = function(object, ntop=700) {
  rv <- matrixStats::rowVars(as.matrix(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
      length(rv)))]
  pca <- prcomp(t(as.matrix(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  names(percentVar) = colnames(pca$x)
  pca$percentVar = percentVar
  return(pca)}

plotFOV = function(eset, metadata, fov_threshold) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  pdat$pcounted = pdat$FovCounted/pdat$FovCount * 100
  ggplot(pdat, aes(rowname, pcounted, color=key2)) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$pcounted))) +
    ylab("percentage of FOV counted") + xlab("Sample.name") +
    geom_hline(yintercept=fov_threshold, color="red") +
    labs(color="pCR")
}
plotFOV.key1 = function(eset, metadata, fov_threshold) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    inner_join(metadata, by=c("rowname"="Sample.name"))
  pdat$pcounted = pdat$FovCounted/pdat$FovCount * 100
  ggplot(pdat, aes(rowname, pcounted, color=key1)) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$pcounted))) +
    ylab("percentage of FOV counted") + xlab("Sample.name") +
    geom_hline(yintercept=fov_threshold, color="red") +
    labs(color="NAC status")
}

plotFOV.BMI = function(eset, metadata, fov_threshold) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  pdat$pcounted = pdat$FovCounted/pdat$FovCount * 100
  ggplot(pdat, aes(rowname, pcounted, color=BMI.Class)) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$pcounted))) +
    ylab("percentage of FOV counted") + xlab("Sample.name") +
    geom_hline(yintercept=fov_threshold, color="red") +
    labs(color="BMI") +
    #scale_colour_gradientn(colours = c("green", "orange", "red"))
    scale_colour_brewer(palette="Reds")
}

plotFOV.invasive = function(eset, metadata, fov_threshold) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  pdat$pcounted = pdat$FovCounted/pdat$FovCount * 100
  ggplot(pdat, aes(rowname, pcounted, color=as.factor(invasive))) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$pcounted))) +
    ylab("percentage of FOV counted") + xlab("Sample.name") +
    geom_hline(yintercept=fov_threshold, color="red") +
    labs(color="Tumour type") +
    #scale_colour_gradient()
    scale_colour_brewer(palette="Set1")
}

plotBD = function(eset, metadata, color_by) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  ggplot(pdat, aes(rowname, BindingDensity, color=key2)) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$BindingDensity))) +
    ylab("Binding density") + xlab("Sample.name") +
    geom_hline(yintercept=0.05, color="red") +
    geom_hline(yintercept=2.25, color="red") +
    labs(color="pCR") +
    scale_colour_brewer(palette="Set1")
}

plotBD.invasive = function(eset, metadata, color_by) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    left_join(metadata, by=c("rowname"="Sample.name"))
  ggplot(pdat, aes(rowname, BindingDensity, color=as.factor(invasive))) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$BindingDensity))) +
    ylab("Binding density") + xlab("Sample.name") +
    geom_hline(yintercept=0.05, color="red") +
    geom_hline(yintercept=2.25, color="red") +
    labs(color="Tumour type") +
    scale_colour_brewer(palette="Set1")
}

plotBD.key1 = function(eset, metadata, color_by) {
  pdat = pData(eset) %>%
    tibble::rownames_to_column() %>%
    inner_join(metadata, by=c("rowname"="Sample.name"))
  ggplot(pdat, aes(rowname, BindingDensity, color=key1)) + geom_point() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x=element_text(size=8)) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = c(0,1.05 * max(pdat$BindingDensity))) +
    ylab("Binding density") + xlab("Sample.name") +
    geom_hline(yintercept=0.05, color="red") +
    geom_hline(yintercept=2.25, color="red") +
    labs(color="NAC status")
}



is_positive = function(column) {
  return(grepl("Pos", column))
}
is_negative = function(column) {
  return(grepl("Neg", column))
}
is_spikein = function(column) {
  return(grepl("Spike", column))
}
is_ligation = function(column) {
  return(grepl("Ligati", column))
}
is_housekeeping = function(column) {
  return(grepl("Housekee", column))
}

extract_pred <- function(eset, predicate, counts=FALSE) {
  if(!counts) {
    counts = data.frame(exprs(eset))
  } else {
    counts = eset
    }
  toplot = counts[predicate(rownames(counts)),] %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("Sample.name", "count", -rowname)
  colnames(toplot) = c("spot", "Sample.name", "count")
  toplot = toplot %>% left_join(metadata, by="Sample.name")
  return(toplot)
}
spotbarplot <- function(toplot) {
  ggplot(toplot,
        aes(Sample.name, count)) + geom_bar(stat='identity') +
    facet_wrap(~spot) +
    theme(axis.text.x = element_blank(),
          text = element_text(size=8))
}
boxplot_expr <- function(eset,predicate) {
    DF <- extract_pred(eset,predicate)%>%
        tidyr::separate(spot,c("a","b","c","d"),sep="_",remove=F)%>%
        tidyr::unite(gene,b,c,sep="_")%>%
        dplyr::select(gene,count)
    ggplot(DF,
        aes(x=gene,y=count)) + geom_boxplot(colour = "black", fill = "#56B4E9")+
        scale_y_continuous(trans = "log2")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
        ylab("counts")
}

hK_mean_exprs_plot <- function(eset, metadata, noise_threshold){
  counts <-  getCounts(eset, metadata)
  hk <-  counts[grepl("Housekeeping", rownames(counts)),]
  hkDF <- as.data.frame(hk)
  tidyHK <- tidyr::gather(hkDF)
  colnames(tidyHK) <- c("sample","count")
  ggplot(tidyHK, aes(sample, count)) + 
    geom_boxplot(colour = "black", fill = "#CCCCCC",outlier.size = 0.5) +
    scale_y_continuous(trans = "log2") + 
    xlab("") + ggtitle("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 5)) +
    geom_smooth(se=T, aes(group=1)) +
    geom_hline(yintercept = (noise_threshold), colour="red")+
    ylab("counts")
}

general_exprs_plot <- function(eset, metadata,noise_threshold){
  counts = getCounts(eset, metadata)
  endG = counts[grepl("Endogenous", rownames(counts)),]
  hk = counts[grepl("Housekeeping", rownames(counts)),]
  all<-rbind(endG,hk)
  allDF<-as.data.frame(all)
  tidyallDF <-tidyr::gather(allDF)
  colnames(tidyallDF)<-c("sample","count")
  ggplot(tidyallDF, aes(sample, count)) + 
    geom_boxplot(colour = "black", fill = "#56B4E9",outlier.size = 0.5) +
    scale_y_continuous(trans = "log2") + 
    xlab("") + ggtitle("") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 5)) +
    geom_smooth(se=T, aes(group=1)) +
    geom_hline(yintercept = noise_threshold,colour="red")
}

