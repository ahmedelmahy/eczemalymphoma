# Define getDEgenes function which takes a DESeqDataSet object
# and returns DESeqDataSet object with significantly expressed genes
getDEgenes <- function(ddsDESeqObject,padj_max = 0.01,lfc_min = 2,
                       verbose = FALSE){
    # change from DESeqDataSet object to DESeqResults object
    res <- results(ddsDESeqObject)
    # Select significantly expressed genes using:
    # 1. log2FoldChange
    #     "The column log2FoldChange is the effect size estimate.
    #     It tells us how much the gene's expression seems to have changed
    #     due to treatment with dexamethasone in comparison to untreated
    #     samples. This value is reported on a logarithmic scale to base 2:
    #     for example, a log2 fold change of 1.5 means that the gene's
    #     expression is increased by a multiplicative factor of 2^1.5≈2.82."
    #
    #     here we use minimum 2 which means only genes that are four times
    #     higher or lower in the treatment group are selected
    #     (log2(4) = 2; 2^2 = 4)
    #
    #2. pvalue
    #     p value indicates the probability that a fold change as strong as the
    #     observed one,or even stronger, would be seen under the situation
    #     described by the null hypothesis.
    res.sel <- res [which(res$padj<padj_max &abs(res$log2FoldChange)>=lfc_min),]
    if (verbose) {
        cat("summary of data when padj_max is ", padj_max)
        print(table(res$padj < padj_max))
    }
    if (verbose) {
        cat("summary of data when lfc_min is ", lfc_min)
        print(table(res$log2FoldChange >= lfc_min))
    }
    if (verbose) {
        cat("summary of data when both are true: ")
        print(table(res$padj<padj_max &abs(res$log2FoldChange)>=lfc_min ))
    }

    # move rownames to listData(just a list of results)  inside
    # res.sel(DESeqResults)
    res.sel$entrez <- rownames(res.sel)
    return(res.sel)
}
#-------------------------------------------------------------------------------
# define getFPKMs which takes:
# 1. ddsDESeqObject(DESeqDataSet) the output of gene expression analysis from
#   DESeq
# 2. ebg (GRangeslist object)
getFPKMs <- function(ddsDESeqObject, ebg, verbose = FALSE) {
    # just make sure colnames the same as ID
    # so now we have:
    #    colnames = 18 # samples
    #    rownames or names = 23057 # genes
    colnames(ddsDESeqObject) <- colData(ddsDESeqObject)$ID
    if (verbose){
        cat("The ddsDESeqObject contains ", length(colnames(ddsDESeqObject)),
            "samples, and ",length(names(ddsDESeqObject)), "genes\n")
    }
    # annotate the gene expression output with information about the genes
    # 1. subset ebg to include only genes in our DESeqDataSet
    ebg.new <- ebg[which(names(ebg)%in%names(ddsDESeqObject)), ]
    # 2. get the rowRanges from the ebg.new
    # so now we have the number of exons and their length for every gene
    rowRanges(ddsDESeqObject) <- ebg.new
    #
    # we use the rowRanges to find the "feature length" which is:
    #    number of basepairs in the union of all GRanges assigned to a given
    #    row of object,
    #    e.g., the union of all basepairs of exons of a given gene.
    # feature length is used to calculate FPKM: fragments per kilobase per
    # million mapped fragments
    mat.norm <- DESeq2::fpkm(ddsDESeqObject, robust=T)
    # mat.norm is a matrix of the  genes x samples but now normalized
    if (verbose) cat("now dimention of fpkm normalized matrix is ",
                     dim(mat.norm))
    return(mat.norm)
}

#-------------------------------------------------------------------------------
# Define getFPKMFCs function to deal with paired data
# which takes:
# 1. mat.norm: a normalized gene expression matrix with fpkm
# 2. ddsDESeqObject: the DESeqDataSet of the same experiment
# returns a matrix of log2(treatment/ normal )

getFPKMFCs <- function(mat.norm, ddsDESeqObject, verbose = FALSE) {
    # remove spaces and use the sample id to be colname for mat.norm
    # [question] here I assume pedigree is just the id of the sample
    colnames(mat.norm) <- gsub("\\s+", "", colData(ddsDESeqObject)$pedigree)
    # split the paired data into ni and i
    mat.ni <- mat.norm[, colData(ddsDESeqObject)$condition=="nI"]
    mat.i <- mat.norm[, colData(ddsDESeqObject)$condition=="I"]
    # make sure mat.ni are the same samples as mat.i
    mat.ni <- mat.ni[, colnames(mat.i)]
    # unify both by dividing
    mat.fcs <- log2((mat.i/mat.ni) + 1)
    if (verbose) cat("\n the input matrix dimentions was ",dim(mat.norm),
                     "\n mat.ni is ",dim(mat.ni),
                     "\n mat.i is ",dim(mat.i), "and the output mat.fcs is",
                     dim(mat.fcs) )
    return(mat.fcs)
}
