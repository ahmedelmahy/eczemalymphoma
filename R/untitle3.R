countData <- matrix(1:100,ncol=4)
condition <- factor(c("A","A","B","B"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)




library(DESeq2)
dat = data
data = dat
dat <- read.csv("df_lym_ecz.csv",header = TRUE ,
                row.names = 1)

class <- factor(dat$class)
dat$class <- NULL

data = t(dat)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = DataFrame(class),
                              design= ~ class)
dds <- DESeq(dds)
res <- results(dds)
res_sel <- res [which(res$padj< config_df$padj_max &
                          abs(res$log2FoldChange)>=config_df$lfc_min),]

genes <- res_sel@rownames
d <- dds[genes,]
d <- as.data.frame(t(counts(d ,
                            normalized = TRUE)))
d$class <- class
return(d)
for (i in 1:10){
    for (j in 1:10){
        print(as.ndata[i,j] < 0)
    }
}
