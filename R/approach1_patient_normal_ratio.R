# a function that takes paired counts matrix


cts <-
cts <- data.frame(t(cts)) # extract counts matrix
cts$sample_ID <- colData(dds.ecz)@listData$ID # get sample ID

cts_ecz <- counts(dds.ecz)
condition_ecz <- paste0("ecz", colData(dds.ecz)@listData$condition)


cts_lym <- counts(dds.lym)
condition_lym <- paste0("lym", colData(dds.lym)@listData$condition)

bind_rows(cts_ecz,cts_lym)
names(cts_lym) <- rownames(cts_lym)
== rownames(cts_ecz)

d <- DGEList(counts=cbind(cts_ecz,cts_lym),group=factor(c(condition_ecz,condition_lym)))


d_lym <- DGEList(counts=cts_lym,group=condition_lym)

d <- edgeR::cbind(d_ecz, d_lym)

# remove less important genes
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

# reset the library size
d$samples$lib.size <- colSums(d$counts)
d$samples

d <- calcNormFactors(d, method = "TMM")
d$samples$norm.factors

#-------------------------------------------------------------------------------




cts_I <- cts[cts$condition == "I",]
cts_nI <- cts[cts$condition == "nI",]


