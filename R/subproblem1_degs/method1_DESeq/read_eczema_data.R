# eczema
# load eczema data
# the data includes:
# 1. dds.test(DESeqDataSet)
# 2. ddsDESeqObject (DESeqDataSet)
# 3. res (data.frame 25673 x 13 )  >> not included in the analysis

load("data/ecz_ddsobj.RData")

# rename ddsDESeqObject
dds.ecz <- ddsDESeqObject

# remove the individual which are both affected skin
dds.ecz <- dds.ecz[,colData(dds.ecz)$pedigree!="10068"]

# train_test = TRUE
# test_now = FALSE
# n_test_ecz = 1
# if(train_test) {
#     dds.ecz <- train_test_split(dds.ecz,n_test_ecz)
# }


# select differential expressed genes
ecz.sel <- getDEgenes(dds.ecz,verbose = T) # lfc and fdr threshhold for degs

# renormalize the data using FPKM (fragments per kilobase per million mapped fragments):
# so that we can compare different experments
ecz.norm <- getFPKMs (dds.ecz, ebg, verbose = TRUE) # fpkm calculations

# unify the paired data into one by creating a new matrix with which equals log2(treatment/ normal )
ecz.fcs <- getFPKMFCs(ecz.norm, dds.ecz,verbose = TRUE)

# table(is.infinite(lym.fcs) | is.na(lym.fcs))
# [question] quality check to change NAs and NAns to 0
ecz.fcs[which(is.infinite(ecz.fcs) | is.na(ecz.fcs))] <- 0
