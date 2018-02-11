cts_ecz_I <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "I"]
d_ecz_I <- DGEList(counts=cts_ecz_I)
d_ecz_I <- calcNormFactors(d_ecz_I, method = "TMM")
v_ecz_I <- as.tibble(t(voom(d_ecz_I)$E))
v_ecz_I$class <- rep("eczI",dim(v_ecz_I)[1])

cts_ecz_nI <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "nI"]
d_ecz_nI <- DGEList(counts=cts_ecz_nI)
d_ecz_nI <- calcNormFactors(d_ecz_nI, method = "TMM")
v_ecz_nI <- as.tibble(t(voom(d_ecz_nI)$E))
v_ecz_nI$class <- rep("ecznI",dim(v_ecz_nI)[1])

cts_lym_I <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "I"]
d_lym_I <- DGEList(counts=cts_lym_I)
d_lym_I <- calcNormFactors(d_lym_I, method = "TMM")
v_lym_I <- as.tibble(t(voom(d_lym_I)$E))
v_lym_I$class <- rep("lymI",dim(v_lym_I)[1])

cts_lym_nI <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "nI"]
d_lym_nI <- DGEList(counts=cts_lym_nI)
d_lym_nI <- calcNormFactors(d_lym_nI, method = "TMM")
v_lym_nI <- as.tibble(t(voom(d_lym_nI)$E))
v_lym_nI$class <- rep("lymnI",dim(v_lym_nI)[1])
#-------------------------------------------------------------------------------
# problem 1
d_train <- rbind(v_ecz_nI, v_lym_nI)
y_train_factor <- as.factor(d_train$class)
d_train$class <- NULL
