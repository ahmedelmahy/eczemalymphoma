load("data/dds_ecz.rda")
load("data/dds_lym.rda")
cts_ecz <- counts(dds_ecz)
cts_lym <- counts(dds_lym)

# combine all to get low expressed genes
d <- cbind(cts_ecz,cts_lym)
d_DGEList <-  DGEList(counts=d)
cpm <- cpm(d_DGEList)
keep <- rowSums(cpm > 3) > 15
keep_genes <- rownames(d_DGEList)[keep]
# using keep_genes
cts_ecz2 <- cts_ecz[keep_genes,]
cts_lym2 <- cts_lym[keep_genes,]
# normalise with voom
v_ecz <- as.data.frame(t(voom(cts_ecz2)$E))
v_lym <- as.data.frame(t(voom(cts_lym2)$E))
# add meta variables
v_ecz$condition <- colData(dds_ecz)@listData$condition
#v_ecz$class <- rep("ecz", dim(cts_ecz)[2])
v_ecz$pedigree <- colData(dds_ecz)@listData$pedigree
v_lym$condition <- colData(dds_lym)@listData$condition
#v_lym$class <- rep("lym", dim(cts_lym)[2])
v_lym$pedigree <- colData(dds_lym)@listData$pedigree

m = v_ecz[v_ecz$condition == "I",]
n = v_ecz[v_ecz$condition == "nI",]
matchh <- match(m$pedigree , n$pedigree)
n = n[matchh,]
sum(n$pedigree != m$pedigree)
m$condition <- NULL
m$pedigree <- NULL
n$condition <- NULL
n$pedigree <- NULL
df_ecz <- m/n
df_ecz$class <- "ecz"
#-------------------------------------------------------------------------------
m = v_lym[v_lym$condition == "I",]
n = v_lym[v_lym$condition == "nI",]
matchh <- match(m$pedigree , n$pedigree)
n = n[matchh,]
sum(n$pedigree != m$pedigree)
m$condition <- NULL
m$pedigree <- NULL
n$condition <- NULL
n$pedigree <- NULL
df_lym <- m/n
df_lym$class <- "lym"
#-------------------------------------------------------------------------------
voom_ratio_df <- as.data.frame(rbind(df_lym, df_ecz))
save(voom_ratio_df,file = "data/voom_ratio_df.rda")
