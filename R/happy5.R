library(cgdsr)
# Create CGDS objec
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
# skcm_tcga is row number 143
getGeneticProfiles(mycgds,"skcm_tcga")[,c(1:2)]
# skcm_tcga_rna_seq_v2_mrna
getCaseLists(mycgds,"skcm_tcga")[,c(1:2)]
# skcm_tcga_rna_seq_v2_mrna
class <- getClinicalData(mycgds,"skcm_tcga_rna_seq_v2_mrna" )$AJCC_PATHOLOGIC_TUMOR_STAGE

df <- getProfileData(mycgds, "NF1", "skcm_tcga_rna_seq_v2_mrna", "skcm_tcga_rna_seq_v2_mrna")
df$class = class
df <- df[class != "",]
library(reshape2)
df_m <- melt(df)
library(ggplot2)
ggplot(df_m, aes(x = class, y = value ))+
    geom_boxplot()
