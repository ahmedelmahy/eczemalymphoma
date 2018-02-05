# load lymphoma data
# the data includes:
# 1. dds.test(DESeqDataSet)
# 2. ddsDESeqObject (DESeqDataSet)
# 3. res (data.frame 25673 x 13 )  >> not included in the analysis
load("data/lym_ddsobj.RData")

# rename ddsDESeqObject
dds.lym <- ddsDESeqObject

# select differential expressed genes
lym.sel <- getDEgenes (dds.lym, verbose = T) # lfc and fdr threshhold for degs



# renormalize the data using FPKM (fragments per kilobase per million mapped fragments):
# so that we can compare different experments
lym.norm <- getFPKMs (dds.lym, ebg, verbose = T) # fpkm calculations

# unify the paired data into one by creating a new matrix with which equals log2(treatment/ normal )
# FCS is for changes
lym.fcs <- getFPKMFCs(lym.norm, dds.lym,verbose = TRUE)
# table(is.infinite(lym.fcs) | is.na(lym.fcs))
# [question] quality check to change NAs and NAns to 0
lym.fcs[which(is.infinite(lym.fcs) | is.na(lym.fcs))] <- 0
