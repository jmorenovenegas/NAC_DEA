## DEA Core script 4. Pathway analysis of the DEA results
# Input: Dataframe with the DEA results obtainedas per coreDEA.R
# Output: table with GSEA analysis.

# Based on `https://github.com/hbc/clark-nanostring`, `https://rawgit.com/hbc/clark-nanostring/master/tumor-set/tumor-set.html#t-test`, `http://bioconductor.org/packages/devel/bioc/vignettes/NanoStringQCPro/inst/doc/vignetteNanoStringQCPro.pdf`

##### Libraries ######
source("0_loadlibraries.R")
orgdb = "org.Hs.eg.db"
biomart_dataset = "hsapiens_gene_ensembl"
keggname = "hsa"
reactomename = "human"
loadpkg("dplyr")
loadpkg("clusterProfiler")
loadpkg("org.Hs.eg.db")
loadpkg("biomaRt")
loadpkg("pathview")
loadpkg("ReactomePA")
loadpkg("DOSE")




mart = biomaRt::useMart(biomart = "ensembl", dataset = biomart_dataset)
entrez = biomaRt::getBM(attributes = c("refseq_mrna", "entrezgene"), mart = mart)
entrez$entrezgene = as.character(entrez$entrezgene)
entrezsymbol = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), mart = mart)
entrezsymbol$entrezgene = as.character(entrezsymbol$entrezgene)

summarize_cp = function(res, comparison) {
    summaries = data.frame()
    for (ont in names(res)) {
        ontsum = summary(res[[ont]])
        ontsum$ont = ont
        summaries = rbind(summaries, ontsum)
    }
    summaries$comparison = comparison
    return(summaries)
}

enrich_cp = function(res, comparison, type="over") {
    res = res %>% data.frame()  %>% left_join(entrez, by = "refseq_mrna") %>% filter(!is.na(entrezgene))
    universe = res$entrezgene 
    if(type=="all"){
        res <- res %>% filter(padj < 0.05)
        genes = res$entrezgene
        # Let the background genes be all the genes not just the genes in the gene expression panel
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
        re = enrichPathway(gene=genes, organism = reactomename, pvalueCutoff=1, qvalueCutoff = 1, pAdjustMethod = "BH")
        # mf = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # cc = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # bp = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # kg = enrichKEGG(gene = genes, universe = universe, organism = keggname, pvalueCutoff = 1, 
        # qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg, re = re)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }
    if(type=="over"){
        res.over <- res %>% filter(log2FoldChange > 1 & padj < 0.05)
        genes = res.over$entrezgene
        # Let the background genes be all the genes not just the genes in the gene expression panel
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
        re = enrichPathway(gene=genes, organism = reactomename, pvalueCutoff=1, qvalueCutoff = 1, pAdjustMethod = "BH")
        # mf = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # cc = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # bp = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # kg = enrichKEGG(gene = genes, universe = universe, organism = keggname, pvalueCutoff = 1, 
        # qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg, re = re)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }

    if(type=="under"){
        res.under <- res %>% filter(log2FoldChange < -1 & padj < 0.05)
        genes = res.under$entrezgene
        mf = enrichGO(genes, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        cc = enrichGO(genes, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        bp = enrichGO(genes, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
        kg = enrichKEGG(gene = genes, organism = keggname, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
        re = enrichPathway(gene=genes, organism = reactomename, pvalueCutoff=1, qvalueCutoff = 1, pAdjustMethod = "BH")
        # mf = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "MF", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # cc = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "CC", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # bp = enrichGO(genes, universe = universe, OrgDb = orgdb, ont = "BP", pAdjustMethod = "BH", 
        # qvalueCutoff = 1, pvalueCutoff = 1)
        # kg = enrichKEGG(gene = genes, universe = universe, organism = keggname, pvalueCutoff = 1, 
        # qvalueCutoff = 1, pAdjustMethod = "BH")
        all = list(mf = mf, cc = cc, bp = bp, kg = kg, re = re)
        all[["summary"]] = summarize_cp(all, comparison)
        return(all)
    }
}

gsea_cp = function(res, comparison) {
    res = res %>% data.frame() %>% left_join(entrez, by = "refseq_mrna") %>% 
        filter(!is.na(entrezgene)) %>% filter(!is.na(log2FoldChange)) %>% filter(!is.na(lfcSE))
    lfc = data.frame(res)[, "log2FoldChange"]
    lfcse = data.frame(res)[, "lfcSE"]
    genes = lfc/lfcse
    names(genes) = res$entrezgene
    genes = genes[order(genes, decreasing = TRUE)]
    cc = gseGO(genes, ont = "CC", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
        pAdjustMethod = "BH", verbose = TRUE)
    mf = gseGO(genes, ont = "MF", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
        pAdjustMethod = "BH", verbose = TRUE)
    bp = gseGO(genes, ont = "bp", OrgDb = orgdb, nPerm = 500, pvalueCutoff = 1, 
        pAdjustMethod = "BH", verbose = TRUE)
    # genes = data.frame(res)[, 'log2FoldChange'] names(genes) = res$entrezgene
    # genes = genes[order(genes, decreasing=TRUE)] genes = genes[!is.na(genes)]
    kg = gseKEGG(geneList = genes, organism = keggname, nPerm = 500, pvalueCutoff = 1, 
        verbose = TRUE)
    re = gsePathway(gene=genes, organism = reactomename, pvalueCutoff=1, pAdjustMethod = "BH")
    if (orgdb == "org.Hs.eg.db") {
        do = DOSE::gseDO(geneList = genes, nPerm = 500, pvalueCutoff = 1, pAdjustMethod = "BH", 
            verbose = TRUE)
        all = list(mf = mf, cc = cc, bp = bp, kg = kg, do = do, re = re)
    } else {
        all = list(mf = mf, cc = cc, bp = bp, kg = kg, re = re)
    }
    all[["summary"]] = summarize_cp(all, comparison)
    return(all)
}


convert_enriched_ids = function(res, entrezsymbol) {
    res = res %>% mutate(geneID = strsplit(as.character(geneID), "/")) %>% tidyr::unnest(geneID) %>% 
        left_join(entrezsymbol, by = c(geneID = "entrezgene")) %>% group_by(ID, 
        Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, ont, 
        comparison) %>% summarise(geneID = paste(geneID, collapse = "/"), symbol = paste(hgnc_symbol, 
        collapse = "/"))
    return(res)
}


convert_core_ids = function(row) {
    core_ids = data.frame(entrezgene = unlist(strsplit(row["core_enrichment"], 
        "/")[[1]])) %>% left_join(entrezsymbol, by = "entrezgene")
    core_symbols = unique(core_ids$hgnc_symbol)
    core_symbols = core_symbols[!is.na(core_symbols)]
    names(core_symbols) = NULL
    return(paste(core_symbols, collapse = "/"))
}

convert_gsea_results = function(res) {
    res$symbols = apply(res, 1, convert_core_ids)
    return(res)
}




