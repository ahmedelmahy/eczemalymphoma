diff.Exp.Genes = c("ESR1", "GATA3", "XBP1", "FOXA1", "ERBB2", "GRB7", "EGFR",
                   "FOXC1", "MYC")

s <- Samples.mRNASeq(format = "csv",
                gene = "ESR1",
                cohort = "SKCM",
                page_size = 2000,
                page = 1)



all.Found = F
page.Counter = 1
mRNA.Exp = list()
page.Size = 2000 # using a bigger page size is faster
while(all.Found == F){
    mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv",
                                               gene = diff.Exp.Genes,
                                               cohort = "BRCA",
                                               tcga_participant_barcode =
                                                   brca.Pats$tcga_participant_barcode,
                                               page_size = page.Size,
                                               page = page.Counter)
    if(nrow(mRNA.Exp[[page.Counter]]) < page.Size)
        all.Found = T
    else
        page.Counter = page.Counter + 1
}
mRNA.Exp = do.call(rbind, mRNA.Exp)
dim(mRNA.Exp)
