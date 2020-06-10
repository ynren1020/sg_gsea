#print("This file is created within Rstudio")

# test data.table::fread function for large data reading ---

library(data.table)
#dataDir <- path.expand("/Users/yren/Projects/sgGSEA")
dataDir <- path.expand("~/Projects/sgGSEA")

# interested cohort ---
cohort <- as.data.frame(data.table::fread(paste0(dataDir, "/PRAD.FPKM.txt")))

# filter out low expressed genes -----
cohort <- cohort[rowMeans(cohort[,-1])>=1, ]
rownames(cohort) <- cohort$GeneID
cohort$GeneID <- NULL
# transpose ---
cohortT <- as.data.frame(t(cohort))
cohortT <- data.frame(patient_id = row.names(cohortT), cohortT)
# interested t_id ---
t_id <- "ENST00000561519"
t_id <- "T280760"
reflnc <- data.table::fread(paste0(dataDir, "/RefLnc_lncRNA_tumor_sample_FPKM_tid"))
mich <- data.table::fread(paste0(dataDir, "/mitranscriptome.expr.fpkm_tid.tsv"))

if (any(grepl(t_id, reflnc[[1]]))) {
dataFiles <- dir(dataDir, pattern = "FPKM$", full.names = TRUE)
meta <- data.table::fread(paste0(dataDir, "/tcga.meta.file.txt"))
# header ---
command0 <- sprintf(
    "grep t_name %s",
    paste(dataFiles, collapse = " ")
)
} else if (any(grepl(t_id, mich[[1]]))) {
    dataFiles <- dir(dataDir, pattern = "fpkm.tsv$", full.names = TRUE) 
    meta <- data.table::fread(paste0(dataDir, "/library_info.txt"))
    # header ---
    command0 <- sprintf(
        "grep transcript_id %s",
        paste(dataFiles, collapse = " ")
    )
} else { message("Can not find the input t_id in the package database.")}


command <- sprintf(
    "grep %s",
    paste(t_id, dataFiles, collapse = " ")
)

# create data.table with header ---
header <- data.table::fread(cmd = command0)
dt <- data.table::fread(cmd = command)
names(dt) <- names(header)

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
    t_dt.join.sub <- t_dt.join[,c("tcga_legacy_sample_id",t_id)]
    tid_cohort <- na.omit(dplyr::full_join(t_dt.join.sub, cohortT,
                                           by = c("tcga_legacy_sample_id" = "patient_id")))
}






