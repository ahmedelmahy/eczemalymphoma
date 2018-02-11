#-------------------------------------------------------------------------------
# I for lesion


cts_ecz_nI <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "nI"]
coldata_ecz_nI <- colData(dds.ecz)[colData(dds.ecz)@listData$condition == "nI",]
coldata_ecz_nI@listData$condition <- rep("ecznI", dim(cts_ecz_nI)[2])

cts_lym_I <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "I"]
coldata_lym_I <- colData(dds.lym)[colData(dds.lym)@listData$condition == "I",]
coldata_lym_I@listData$condition <- rep("lymI", dim(cts_lym_I)[2])

cts_lym_nI <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "nI"]
coldata_lym_nI <- colData(dds.lym)[colData(dds.lym)@listData$condition == "nI",]
coldata_lym_nI@listData$condition <- rep("lymnI", dim(cts_lym_nI)[2])


# cts.ecz <- cts.ecz[apply(cts.ecz, 1, function(x) sum(x == 0)) < (ncol(cts.ecz)* 0.8),]
#x <- DGEList(cts.ecz)
condition=factor(cts.ecz$condition)
design=model.matrix(~0+condition)

v=voom(cts.ecz)$E
cor=duplicateCorrelation(v,design,block=info$individual)

cts.ecz$class <- rep("lym", dim(cts.ecz)[1])

load("data/lym_ddsobj.RData") # load data object
dds.lym <- ddsDESeqObject # rename ddsDESeqObject
cts.lym <- as.tibble(t(counts(dds.lym, normalized = FALSE))) # extract counts matrix
cts.lym$ID <- colData(dds.lym)@listData$ID # get sample ID
cts.lym$condition <- colData(dds.lym)@listData$condition
cts.lym$class <- rep("lym", dim(cts.lym)[1])

# join both



# clean env
rm(dds.test,res,dds.ecz,dds.lym, ddsDESeqObject,cts.lym,cts.ecz)
#-------------------------------------------------------------------------------
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
#removal of samples with expression estimates with counts in less than 20% of cases.





