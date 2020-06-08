

sg_gsea <- function(tid_cohort, method = "logFC", geneset, pathway){
    # tid_cohort <- reflncT.tid.pradT
    
    if (method = "logFC") {
    # by group, calculate mean for each gene ---
    res <- aggregate(. ~ group, data = tid_cohort[, -c(1:2)], FUN = mean)
    rownames(res) <- res$group
    res$group <- NULL
    resT <- as.data.frame(t(res))
    
    # log2 FC of high to low for each gene ---
    resT$FC <- resT$high/resT$low
    resT$log2FC <- log2(resT$FC)
    
    # rownames to genes
    resT$genes <- rownames(resT)
    resT <- resT[ ,c(5,4)]
    resT <-resT[order(-resT[,2]),] # log2FC decrease 
    resT
    
    # gsea using fgsea -----
    ranks <- tibble::deframe(resT)
    #head(ranks, 20)
    } else {
        # correlation -----
        res <- NULL
        res <- cor(tid_cohort[,1], tid_cohort[,3:ncol(tid_cohort)],
                   use = 'pairwise.complete.obs')
        
        # make res ready for gsea -----
        res <- as.data.frame(t(as.data.frame(res)))
        names(res) <- tid
        res$genes <- rownames(res)
        res <- res[ ,c(2,1)]
        res<-res[order(-res[,2]),] # cor decrease 
        res
        
        # gsea using fgsea -----
        ranks <- tibble::deframe(res)
    }
    # Load the pathways into a named list-----
    pathways.hallmark <- fgsea::gmtPathways(geneset) # gmt 
    
    fgseaRes <- fgsea::fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
    
     # tidy the results -----
    fgseaResTidy <- fgseaRes[order(-NES,padj),]
    print(fgseaResTidy)
    data.table::fwrite(fgseaResTidy, paste0(tid, "_", strsplit(cohort,"[.]")[[1]][1], ".txt"),
                       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    
    # gsea plot -----
    #args5 <- "HALLMARK_ANDROGEN_RESPONSE"
    #args5 <- "HALLMARK_MYC_TARGETS_V1"
    #args5 <- NULL
    if(!is.null(pathway)){
        pdf(file = paste0(pathway, "_", tid, ".pdf"))
        fgsea::plotEnrichment(pathways.hallmark[[pathway]], ranks) + ggplot2::labs(title = pathway)}
    if(!is.null(pathway)){
        dev.off()
    }
    
}