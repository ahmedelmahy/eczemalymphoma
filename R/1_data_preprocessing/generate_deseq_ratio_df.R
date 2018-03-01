
#dds_ecz <- dds_ecz[which(colData(dds_ecz)@listData$pedigree != "10068"),]

library(DESeq2)
load("data/dds_ecz.rda")
load("data/dds_lym.rda")
load("data/ebg.rda")
source("R/1_data_preprocessing/DESeq_processing_functions.R")


ecz_sel <- getDEgenes(dds_ecz,verbose = T)
ecz_norm <- getFPKMs (dds_ecz, ebg, verbose = TRUE)
ecz_fcs <- getFPKMFCs(ecz_norm, dds_ecz,verbose = TRUE)
ecz_fcs[which(is.infinite(ecz_fcs) | is.na(ecz_fcs))] <- 0


lym_sel <- getDEgenes (dds_lym, verbose = T)
lym_norm <- getFPKMs (dds_lym, ebg, verbose = T)
lym_fcs <- getFPKMFCs(lym_norm, dds_lym,verbose = TRUE)
lym_fcs[which(is.infinite(lym_fcs) | is.na(lym_fcs))] <- 0

gene_mat <- merge(data.frame(lym_sel), data.frame(ecz_sel), by="entrez", all=T)
lym_mat <- data.frame(lym_fcs[which(rownames(lym_fcs)%in%gene_mat$entrez), ])
ecz_mat <- data.frame(ecz_fcs[which(rownames(ecz_fcs)%in%gene_mat$entrez), ])
data_mat <- cbind(lym_mat, ecz_mat[rownames(lym_mat), ])
data_mat <- t(data_mat[which(apply(data_mat, 1,
                                   function(x) length(which(
                                       is.infinite(x) | is.na(x))))==0), ])

class_labels <- as.factor(c(rep("lym", ncol(lym_mat)),
                            rep("ecz", ncol(ecz_mat))))
deseq_ratio_df <- data.frame(data_mat)
deseq_ratio_df$class <- class_labels



save(deseq_ratio_df,file = "data/deseq_ratio_df.rda")


