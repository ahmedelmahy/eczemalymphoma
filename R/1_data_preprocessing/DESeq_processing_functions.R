# Define getDEgenes function which takes a DESeqDataSet object
# and returns DESeqDataSet object with significantly expressed genes
getDEgenes <- function(ddsDESeqObject,padj_max = 0.01,lfc_min = 2,
                       verbose = FALSE){
    # change from DESeqDataSet object to DESeqResults object
    res <- results(ddsDESeqObject)
    res_sel <- res [which(res$padj<padj_max &abs(res$log2FoldChange)>=lfc_min),]
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
    # res_sel(DESeqResults)
    res_sel$entrez <- rownames(res_sel)
    return(res_sel)
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
    ebg_new <- ebg[which(names(ebg)%in%names(ddsDESeqObject)), ]
    # 2. get the rowRanges from the ebg_new
    # so now we have the number of exons and their length for every gene
    rowRanges(ddsDESeqObject) <- ebg_new
    #
    # we use the rowRanges to find the "feature length" which is:
    #    number of basepairs in the union of all GRanges assigned to a given
    #    row of object,
    #    e.g., the union of all basepairs of exons of a given gene.
    # feature length is used to calculate FPKM: fragments per kilobase per
    # million mapped fragments
    mat_norm <- DESeq2::fpkm(ddsDESeqObject, robust=T)
    # mat_norm is a matrix of the  genes x samples but now normalized
    if (verbose) cat("now dimention of fpkm normalized matrix is ",
                     dim(mat_norm))
    return(mat_norm)
}


getFPKMFCs <- function(mat_norm, ddsDESeqObject, verbose = FALSE) {
    colnames(mat_norm) <- gsub("\\s+", "", colData(ddsDESeqObject)$pedigree)
    mat_ni <- mat_norm[, colData(ddsDESeqObject)$condition=="nI"]
    mat_i <- mat_norm[, colData(ddsDESeqObject)$condition=="I"]
    # make sure mat_ni are the same samples as mat_i
    mat_ni <- mat_ni[, colnames(mat_i)]
    # unify both by dividing
    #mat_fcs <- log2((mat_i/mat_ni) + 1)
    mat_fcs <- (scale(mat_i)-scale(mat_ni))   # worked not bad
    #mat_fcs <- (scale(mat_i)/scale(mat_ni))   # worked bad

    if (verbose) cat("\n the input matrix dimentions was ",dim(mat_norm),
                     "\n mat_ni is ",dim(mat_ni),
                     "\n mat_i is ",dim(mat_i), "and the output mat_fcs is",
                     dim(mat_fcs) )
    return(mat_fcs)
}
#
# #dds <- dds_ecz
# #n_test = n_test_ecz
train_test_split <- function(dds, n_test){
    if(test_now) {
        # reproduce the DESeq results
        # get the counts
        cts <- counts(dds)
        coldata <- colData(dds)
        if (!exists("test_s")) stop("you should train first")
        # include the test part
        coldata <- coldata[test_s,]
        cts <- cts[,test_s]

    }else{
        # reproduce the DESeq results
        # get the counts
        cts <- counts(dds, normalized = FALSE)
        coldata <- colData(dds)
        # define test_s as a global variable
        test_s <<- which(coldata$pedigree %in% sample(unique(coldata$pedigree),
                                                      n_test))
        # exclude the test part
        coldata <- coldata[-test_s,]
        cts <- cts[,-test_s]
    }
    dds0 <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design= ~ pedigree + condition )
    library("BiocParallel")
    dds0 <- DESeq(dds0,full = ~ pedigree + condition, reduced = ~pedigree,
                  test="LRT",parallel = TRUE, BPPARAM=MulticoreParam(24))
    return(dds0)
}





