#' Single lncRNA Gene Set Enrichment Analysis
#' 
#' TCGA samples are either stratified by lncRNA's expression and then calculate
#' logFC of high vs low express groups, or calculate correlation of genes' expression
#' with lncRNA's expression, pre-ranked file will be used for pre ranked GSEA using fgsea
#' 
#' @param tid_cohort A output of data frame from pre_gsea function
#' @param t_id A string of character, showing transcript id, same as pre_gsea
#' @param method A string of character, showing which metric is used for rank
#' @param cohort A string of character, showing one of TCGA studies
#' @param geneset A string of character, showing the path where gmt is 
#' @param pathway A string of character, showing one of pathway in the gmt file
#' 
#' @import tibble
#' @import fgsea
#' @import data.table
#' 
#' @return A dataframe showing the GSEA results,if pathway is not null, its enrichment plot 
#'         will be produced too.
#'         
#' @example 
#' sg_gsea(tid_cohort = test, t_id = "ENST00000430998", cohort = "BRCA", method = "logFC",
#'         geneset = "/Users/yren/Projects/sgGSEA/gmt/h.all.v7.0.symbols.gmt",
#'         pathway = NULL)
#' @export
sg_gsea <- function(cohort, t_id, method, geneset, pathway = NULL)
    {
    tid_cohort <- pre_gsea(cohort, t_id)
    
    if (method == "logFC") {
    tid_cohort$group <- ifelse(tid_cohort[[2]] > median(tid_cohort[[2]]), "high", "low")
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
        res <- cor(tid_cohort[,2], tid_cohort[,3:ncol(tid_cohort)],
                   use = 'pairwise.complete.obs')
        
        # make res ready for gsea -----
        res <- as.data.frame(t(as.data.frame(res)))
        names(res) <- t_id
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
    data.table::fwrite(fgseaResTidy, paste0(t_id, "_", cohort, ".txt"),
                       col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    
    # gsea plot -----
    #pathway <- "HALLMARK_ANDROGEN_RESPONSE"
    #pathway <- "HALLMARK_MYC_TARGETS_V1"
    #pathway <- NULL
    if(!is.null(pathway)){
        pdf(file = paste0(pathway, "_", t_id, ".pdf"))
        enrichplot <- fgsea::plotEnrichment(pathways.hallmark[[pathway]], ranks) + ggplot2::labs(title = pathway)
        print(enrichplot)
        }
    if(!is.null(pathway)){
        
        dev.off()
   }
}

