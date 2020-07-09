#' process expression data for sg_GSEA function
#' 
#' Find TCGA samples appear both in RNAseq expression (FPKM) data frame and lncRNAs 
#' data frame, if the sample ids are different, use tcga.meta.info to match. 
#' 
#' @param cohort A character string which is the short name for cohorts in TCGA
#' @param t_id A character string of transcript id, either start with ENST or T
#' 
#' @return A data frame of lncRNA's transcript and other genes' expression in interested cohort
#' 
#' @import data.table
#' @import stringr
#' @import dplyr
#' 
#' @details cohort should be one of studies in TCGA, for example, "PRAD" can be used as 
#'          the function's input to stand for prostate cancer; t_id can be from mitranscriptome
#'          database or reflnc database. 
#' 
#' @example 
#' pre_gsea("PRAD", "T280760")
#' pre_gsea("BRCA","ENST00000430998")
#' 
#' @export
pre_gsea <- function(cohort, t_id){
 
    dataDir <- path.expand("~/Projects/sgGSEA")
    
    # interested cohort ---
    cohort <- as.data.frame(data.table::fread(paste0(dataDir, 
                                                     "/", cohort, ".FPKM.txt")))
    
    # filter out low expressed genes -----
    cohort <- cohort[rowMeans(cohort[,-1])>=1, ]
    rownames(cohort) <- cohort$GeneID
    cohort$GeneID <- NULL
    # transpose ---
    cohortT <- as.data.frame(t(cohort))
    cohortT <- data.frame(patient_id = row.names(cohortT), cohortT)
    # interested t_id ---
    #t_id <- "ENST00000561519"
    #t_id <- "T280760"
    reflnc <- data.table::fread(paste0(dataDir, "/RefLnc_lncRNA_tumor_sample_FPKM_tid"))
    mich <- data.table::fread(paste0(dataDir, "/mitranscriptome.expr.fpkm_tid.tsv"))
    
    if (any(grepl(t_id, reflnc[[1]]))) {
        dataFiles <- dir(dataDir, pattern = "FPKM$", full.names = TRUE)
        meta <- data.table::fread(paste0(dataDir, "/tcga.meta.file.txt"))
        
    } else if (any(grepl(t_id, mich[[1]]))) {
        dataFiles <- dir(dataDir, pattern = "fpkm.tsv.gz$", full.names = TRUE) 
        meta <- data.table::fread(paste0(dataDir, "/library_info.txt"))
       
    } else { message("Can not find the input t_id in the package database.
                     Users should provide their lncRNAs expression matrix for 
                     TCGA samples.")}
    
    # create data.table with header ---
    dt <- data.table::fread(dataFiles)
    names(dt)[1] <- "transcript_id"
    dt <- dt[dt$transcript_id==t_id, ]
    
    # transpose ---
    dt <- as.data.frame(dt)
    rownames(dt) <- dt[,1]
    dt[,1] <- NULL
    t_dt <- as.data.frame(t(dt))
    
    # test if sample name end with "gdc_realn_rehead"
    if (any(grepl("gdc_realn_rehead", rownames(t_dt)))){
        t_dt$library_id <- stringr::str_split(rownames(t_dt),"_",simplify = TRUE)[,1]
        t_dt.join <- dplyr::left_join(t_dt, meta, by = "library_id")
        t_dt.join.sub <- t_dt.join[,c(3,1)]
        tid_cohort <- na.omit(dplyr::full_join(t_dt.join.sub, cohortT,
                                               by = c("ALIQUOT_BARCODE" = "patient_id")))
    } else {
        t_dt$library_id <- rownames(t_dt)
        t_dt.join <- dplyr::left_join(t_dt, meta, by = "library_id")
        #t_dt.join.sub <- t_dt.join[,c("tcga_legacy_sample_id",t_id)]
        t_dt.join.sub <- t_dt.join[,c(13,1)]
        tid_cohort <- na.omit(dplyr::full_join(t_dt.join.sub, cohortT,
                                               by = c("tcga_legacy_sample_id" = "patient_id")))
    }
    
    names(tid_cohort)[1] <- "patient_id"
    return(tid_cohort)
    
        
}