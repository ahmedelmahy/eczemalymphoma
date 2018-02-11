# I for lesion
cts_ecz_I <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "I"]
coldata_ecz_I <- colData(dds.ecz)[colData(dds.ecz)@listData$condition == "I",]
coldata_ecz_I@listData$condition <- rep("eczI", dim(cts_ecz_I)[2])

cts_ecz_nI <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "nI"]
coldata_ecz_nI <- colData(dds.ecz)[colData(dds.ecz)@listData$condition == "nI",]
coldata_ecz_nI@listData$condition <- rep("ecznI", dim(cts_ecz_nI)[2])

cts_lym_I <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "I"]
coldata_lym_I <- colData(dds.lym)[colData(dds.lym)@listData$condition == "I",]
coldata_lym_I@listData$condition <- rep("lymI", dim(cts_lym_I)[2])

cts_lym_nI <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "nI"]
coldata_lym_nI <- colData(dds.lym)[colData(dds.lym)@listData$condition == "nI",]
coldata_lym_nI@listData$condition <- rep("lymnI", dim(cts_lym_nI)[2])

#-------------------------------------------------------------------------------
cts_eczI_lymI <- cbind(cts_ecz_I,cts_lym_I)
coldata_eczI_lymI <- rbind(coldata_ecz_I,coldata_lym_I)


cts_eczI_lymnI <- cbind(cts_ecz_I,cts_lym_nI)
coldata_eczI_lymnI <- rbind(coldata_ecz_I,coldata_lym_nI)


cts_ecznI_lymI <- cbind(cts_ecz_nI,cts_lym_I)
coldata_ecznI_lymI <- rbind(coldata_ecz_nI,coldata_lym_I)

cts_ecznI_lymnI <- cbind(cts_ecz_nI,cts_lym_nI)
coldata_ecznI_lymnI <- rbind(coldata_ecz_nI,coldata_lym_nI)

#-------------------------------------------------------------------------------
dds_eczI_lymI <- DESeqDataSetFromMatrix(countData = cts_eczI_lymI,
                               colData = coldata_eczI_lymI,
                               design= ~ condition )

dds_eczI_lymnI <- DESeqDataSetFromMatrix(countData = cts_eczI_lymnI,
                               colData = coldata_eczI_lymnI,
                               design= ~ condition )

dds_ecznI_lymI <- DESeqDataSetFromMatrix(countData = cts_ecznI_lymI,
                               colData = coldata_ecznI_lymI,
                               design= ~ condition )

dds_ecznI_lymnI <- DESeqDataSetFromMatrix(countData = cts_ecznI_lymnI,
                               colData = coldata_ecznI_lymnI,
                               design= ~ condition )

#-------------------------------------------------------------------------------
dds_eczI_lymI <- DESeq(dds_eczI_lymI)
dds_eczI_lymnI <- DESeq(dds_eczI_lymnI)
dds_ecznI_lymI <- DESeq(dds_ecznI_lymI)
dds_ecznI_lymnI<- DESeq(dds_ecznI_lymnI)

#-------------------------------------------------------------------------------

dds_eczI_lymI_deg <- getDEgenes(dds_eczI_lymI, padj_max = 0.01,lfc_min = 3,
                                verbose = TRUE)
dds_eczI_lymnI_deg <- getDEgenes(dds_eczI_lymnI,padj_max = 0.01,lfc_min = 2,
                                 verbose = TRUE)
dds_ecznI_lymI_deg <- getDEgenes(dds_ecznI_lymI,padj_max = 0.01,lfc_min = 2,
                                 verbose = TRUE)
dds_ecznI_lymnI_deg <- getDEgenes(dds_ecznI_lymnI, padj_max = 0.01,lfc_min = 1,
                                  verbose = TRUE)


eIlI <- dds_eczI_lymI_deg$entrez
eIlnI <- dds_eczI_lymnI_deg$entrez
enIlI <- dds_ecznI_lymI_deg$entrez
enIlnI <- dds_ecznI_lymnI_deg$entrez

message("variables eIlI, eIlnI, enIlI, enIlnI")
#i <- setdiff(eIlI,eIlnI)

#MDS_for_selected_genes(dds.ecz = dds.ecz,dds.lym = dds.lym,sel_genes =i)

#-------------------------------------------------------------------------------
#
# dds_eczI_lymI_norm <- getFPKMs(dds_eczI_lymI,ebg, verbose = TRUE)
# dds_eczI_lymnI_norm <- getFPKMs(dds_eczI_lymnI,ebg, verbose = TRUE)
# dds_ecznI_lymI_norm <- getFPKMs(dds_ecznI_lymI,ebg, verbose = TRUE)
# dds_ecznI_lymnI_norm <- getFPKMs(dds_ecznI_lymnI,ebg, verbose = TRUE)
#
# #------------------------------------------------------------------------------

