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
