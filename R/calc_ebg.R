# hg19GeneLengths <- function(symbols)
# {
#     require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#     require(org.Hs.eg.db)
#     exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
#     egs    = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
#     sapply(egs,function(eg)
#     {
#         exons = exons.db[[eg]]
#         exons = reduce(exons)
#         sum( width(exons) )
#     })
# }
# hg19GeneLengths( c('STAT1','CXCL10','ACTB','PDCD1') )

#
x = i[1]
ebg['8631']
# library("GenomicFeatures")
#
# txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.91.gtf", format = "gtf", circ_seqs = character())
# txdb
# ebg <- exonsBy(txdb, by="gene")
# ebg
