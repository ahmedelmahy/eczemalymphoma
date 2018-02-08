#-------------------------------------------------------------------------------
# Get the counts matrix from DESeqDataSet
cts <- counts(dds.lym, normalized = FALSE)
# # # below replicates Richa's DESeq results to make sure cts are correct
# # coldata <- colData(dds.lym)
# # dds0 <- DESeqDataSetFromMatrix(countData = cts,
# #                               colData = coldata,
# #                               design= ~ pedigree + condition )
# dds0 <- DESeq(dds0,full = ~ pedigree + condition, reduced = ~pedigree,
#                           test="LRT", parallel=TRUE, BPPARAM=MulticoreParam(3))
# # lym.sel0 <- getDEgenes (dds0, verbose = T) # lfc and fdr threshhold for degs
# # lym.norm0 <- getFPKMs (dds0, ebg, verbose = T) # fpkm calculations
# # lym.fcs0 <- getFPKMFCs(lym.norm0, dds0,verbose = TRUE)
# # lym.fcs0[which(is.infinite(lym.fcs0) | is.na(lym.fcs0))] <- 0
# # sum(lym.fcs0 != lym.fcs) == 0
#-------------------------------------------------------------------------------

