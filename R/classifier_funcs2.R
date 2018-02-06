getNormCounts <- function(ddsDESeqObject) {
  dds.i <- ddsDESeqObject[, colData(ddsDESeqObject)$condition=="I"]
  dds.ni <- ddsDESeqObject[, colData(ddsDESeqObject)$condition=="nI"]
  sf.ni <- estimateSizeFactors(dds.ni)
  sf.i <- estimateSizeFactors(dds.i)
  norm.counts <- counts(dds.i, normalized=T) / counts(dds.ni, normalized=T)
  return(norm.counts)
}
tpm = function (counts, effective_lengths) {
  rate = log(counts) - log(effective_lengths)
  exp(rate - log(sum(exp(rate))) + log(10 ^ 6))
}
getFPKMs <- function(ddsDESeqObject, ebg) {
  colnames(ddsDESeqObject) <- colData(ddsDESeqObject)$ID
  ebg.new <- ebg[which(names(ebg)%in%names(ddsDESeqObject)==T), ]
  rowRanges(ddsDESeqObject) <- ebg.new
  mat.norm <- fpkm(ddsDESeqObject, robust=T)  
  return(mat.norm)
}
getDEgenes <- function(ddsDESeqObject){
  res <- results(ddsDESeqObject)
  res.sel <- res [which(res$padj<0.01 &abs(res$log2FoldChange)>=2), ]
  res.sel$entrez <- rownames(res.sel)
  return(res.sel)
}
getFPKMFCs <- function(mat.norm, ddsDESeqObject) {
  colnames(mat.norm) <- gsub("\\s+", "", colData(ddsDESeqObject)$pedigree)
  mat.ni <- mat.norm[, colData(ddsDESeqObject)$condition=="nI"]
  mat.i <- mat.norm[, colData(ddsDESeqObject)$condition=="I"]
  mat.ni <- mat.ni[, colnames(mat.i)]
  mat.fcs <- log2((mat.i/mat.ni) + 1)
  return(mat.fcs)
}

plotFeatures <- function(mat.norm, dds, features, out.file){
  pdf(out.file)
  my.features <- data.frame(cbind(t(mat.norm[which(rownames(mat.norm)%in%features), ]),
                                  as.matrix(colData(dds)$condition), as.matrix(colData(dds)$ID)))
  names(my.features) <- c(features, "condition", "ID")
  plot.mat <- melt(features, id=c("condition", "ID"))
  plot.mat$value <- as.numeric(plot.mat$value)
  ggplot(plot.mat, aes(y=value,x=variable, fill=condition)) +geom_boxplot()
  bar.cols <- c("gray", "white")
  lapply(1:length(gene.features), FUN=function(i){
    barplot(plot.mat$value[which(plot.mat$variable==features[i])], col =bar.cols[plot.mat$condition], main=features[i])
  })
  dev.off()
}