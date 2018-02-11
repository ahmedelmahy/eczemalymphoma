rm(new.data)
sel_genes <- i
MDS_for_selected_genes <- function(dds.ecz, dds.lym, sel_genes_list,titles,
                                   include_only_disease = FALSE,
                                   measurement = "cpm",ebg = NA,...){
    if (length(sel_genes_list) != length(titles)) {stop("titles are not the same
                                                       length as selected_genes_list")}
    l <- length(titles)
    if( l== 1)par(mfrow = c(1,1))
    if( l== 2)par(mfrow = c(1,2))
    if( l %in% 3:4)par(mfrow = c(2,2))
    if( l %in% 5:6)par(mfrow = c(2,3))
    if( l %in% 7:9)par(mfrow = c(3,3))
    if( l %in% 10:13)par(mfrow = c(3,4))

    i = 1
    for(sel_genes in sel_genes_list){
        if(include_only_disease){
            cts_ecz <- counts(dds.ecz)[rownames(counts(dds.ecz)) %in% sel_genes,
                                       colData(dds.ecz)@listData$condition == "I"]
            coldata_ecz_I <- colData(dds.ecz)[colData(dds.ecz)@listData$condition == "I",]
            coldata_ecz_I@listData$condition <- rep("eczI", dim(cts_ecz_I)[2])

            cts_lym <- counts(dds.lym)[rownames(counts(dds.lym)) %in% sel_genes,
                                       colData(dds.lym)@listData$condition == "I"]
            coldata_lym_I <- colData(dds.lym)[colData(dds.lym)@listData$condition == "I",]
            coldata_lym_I@listData$condition <- rep("lymI", dim(cts_lym_I)[2])

            group=factor(c(coldata_ecz_I@listData$condition,
                           coldata_lym_I@listData$condition))
            d <- DGEList(counts=cbind(cts_ecz,cts_lym),
                         group = group)
        }else{
            cts_ecz <- counts(dds.ecz)[rownames(counts(dds.ecz)) %in% sel_genes,]
            condition_ecz <- paste0("ecz", colData(dds.ecz)@listData$condition)

            cts_lym <- counts(dds.lym)[rownames(counts(dds.lym)) %in% sel_genes,]
            condition_lym <- paste0("lym", colData(dds.lym)@listData$condition)

            group=factor(c(condition_ecz,condition_lym))
            d <- DGEList(counts=cbind(cts_ecz,cts_lym),
                         group = group)
        }
        d <- calcNormFactors(d, method = "TMM")
        #d$samples$norm.factors
        # MDS plot

        col.group <- group
        levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
        col.group <- as.character(col.group)

        if (measurement == "rpkm"){
            if (length(ebg) < 1) stop("add ebg parameter")
            # ebg.new <- ebg[which(names(ebg)%in%most_genes), ]
            # sapply(1:length(ebg.new), function(x) {cat(x);
            #     ebg[[x]]@ranges@width})
            gene_length <- sapply(sel_genes,function(x) sum(ranges(ebg[x])[[1]]@width))
            rpkm_selected_genes <- rpkm(d, gene.length = gene_length)
            plotMDS(rpkm_selected_genes, labels=group, col=col.group,...)
        } else {
            lcpm <- cpm(d, log=TRUE)
            plotMDS(lcpm, labels=group, col=col.group,...)
        }


        title(titles[i])
        i <- i + 1
    }
    }
