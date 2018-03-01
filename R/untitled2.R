library(magrittr)

load("data/dds_ecz.rda")
load("data/dds_lym.rda")

# generate a dataframe for shiny app
cts_ecz_I <- counts(dds_ecz)[,colData(dds_ecz)@listData$condition == "I"]
cts_lym_I <- counts(dds_lym)[,colData(dds_lym)@listData$condition == "I"]

# assert genes
dim(cts_ecz_I) == dim(cts_lym_I)
sum(rownames(cts_ecz_I) != rownames(cts_lym_I))

#
df <- cbind(cts_ecz_I, cts_lym_I)
df_final <- as.data.frame(t(df))
df_final %>% dim
df_final$class <- c(rep("ecz",dim(cts_ecz_I)[2]),
                    rep("lym",dim(cts_lym_I)[2]))

colnames(df_final)
table(sapply(df_final,class))
#
write.csv(df_final, file = "df_lym_ecz.csv")
