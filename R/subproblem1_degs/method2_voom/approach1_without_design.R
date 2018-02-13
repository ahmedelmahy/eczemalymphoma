cts_ecz <- counts(dds.ecz)
cts_lym <- counts(dds.lym)

# combine all to get low expressed genes
d <- cbind(cts_ecz,cts_lym)
d_DGEList <-  DGEList(counts=d)
cpm <- cpm(d_DGEList)
keep <- rowSums(cpm > .5) > 3
table(keep)

keep_genes <- rownames(d_DGEList)[keep]

# using keep_genes
cts_ecz2 <- cts_ecz[keep_genes,]
cts_lym2 <- cts_lym[keep_genes,]
# normalise with voom
v_ecz <- as.data.frame(t(voom(cts_ecz2)$E))
v_lym <- as.data.frame(t(voom(cts_lym2)$E))

# add meta variables
v_ecz$condition <- colData(dds.ecz)@listData$condition
#v_ecz$class <- rep("ecz", dim(cts_ecz)[2])
v_ecz$pedigree <- colData(dds.ecz)@listData$pedigree

v_lym$condition <- colData(dds.lym)@listData$condition
#v_lym$class <- rep("lym", dim(cts_lym)[2])
v_lym$pedigree <- colData(dds.lym)@listData$pedigree


m = v_ecz[v_ecz$condition == "I",]
n = v_ecz[v_ecz$condition == "nI",]
matchh <- match(m$pedigree , n$pedigree)
n = n[matchh,]
sum(n$pedigree != m$pedigree)
m$condition <- NULL
m$pedigree <- NULL
n$condition <- NULL
n$pedigree <- NULL
df_ecz <- m/n
df_ecz$class <- "ecz"
#-------------------------------------------------------------------------------
m = v_lym[v_lym$condition == "I",]
n = v_lym[v_lym$condition == "nI",]
matchh <- match(m$pedigree , n$pedigree)
n = n[matchh,]
sum(n$pedigree != m$pedigree)
m$condition <- NULL
m$pedigree <- NULL
n$condition <- NULL
n$pedigree <- NULL
df_lym <- m/n
df_lym$class <- "lym"
#-------------------------------------------------------------------------------
mat <- as.data.frame(rbind(df_lym, df_ecz))


# Now we will change df to be new.data
# and try caret again


table(v_lym$lym_pedigree)
cts_ecz_I <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "I"]
d_ecz_I <- DGEList(counts=cts_ecz_I)
d_ecz_I <- calcNormFactors(d_ecz_I, method = "TMM")
v_ecz_I <- as.data.frame(t(voom(d_ecz_I)$E))
v_ecz_I$class <- rep("eczI",dim(v_ecz_I)[1])

cts_ecz_nI <- counts(dds.ecz)[,colData(dds.ecz)@listData$condition == "nI"]
d_ecz_nI <- DGEList(counts=cts_ecz_nI)
d_ecz_nI <- calcNormFactors(d_ecz_nI, method = "TMM")
v_ecz_nI <- as.data.frame(t(voom(d_ecz_nI)$E))
v_ecz_nI$class <- rep("ecznI",dim(v_ecz_nI)[1])

cts_lym_I <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "I"]
d_lym_I <- DGEList(counts=cts_lym_I)
d_lym_I <- calcNormFactors(d_lym_I, method = "TMM")
v_lym_I <- as.data.frame(t(voom(d_lym_I)$E))
v_lym_I$class <- rep("lymI",dim(v_lym_I)[1])

cts_lym_nI <- counts(dds.lym)[,colData(dds.lym)@listData$condition == "nI"]
d_lym_nI <- DGEList(counts=cts_lym_nI)
d_lym_nI <- calcNormFactors(d_lym_nI, method = "TMM")
v_lym_nI <- as.data.frame(t(voom(d_lym_nI)$E))
v_lym_nI$class <- rep("lymnI",dim(v_lym_nI)[1])
#-------------------------------------------------------------------------------
# problem 1
# predict the person
# nI nI
d_train <- rbind(v_ecz_nI, v_lym_nI)
y_train_factor <- as.factor(d_train$class)
d_train$class <- NULL


RF_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                        y_train_factor = y_train_factor,
                                                        model_method = "rf",
                                                        fitControl = fitControl,
                                                        sampling_method = "smote",verboseIter = TRUE)


RF_smote_sampling_selected_genes <- lapply(1:10, FUN=function(i) {
    print(RF_smote_sampling$bestTune$mtry)
    RF_smote_sampling <-
        balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                           y_train_factor = y_all_data,
                                           model_method = "rf",
                                           fitControl = fitControl,
                                           sampling_method = "smote",
                                           tuneGrid =
                                               data.frame(mtry=RF_smote_sampling$bestTune$mtry),
                                           verboseIter = TRUE)
    var.imp.rf <- varImp(RF_smote_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance),
                                   var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)),
                                   decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
})

rf_smote_sampling_selected_genes <- sort(table(unlist(lapply(
    RF_smote_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]
RF_smote_selected_genes_nI_nI <- gsub("X","",names(rf_smote_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------
# nI nI
d_train <- rbind(v_ecz_I, v_lym_nI)
y_train_factor <- as.factor(d_train$class)
d_train$class <- NULL


RF_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                        y_train_factor = y_train_factor,
                                                        model_method = "rf",
                                                        fitControl = fitControl,
                                                        sampling_method = "smote",verboseIter = TRUE)


RF_smote_sampling_selected_genes <- lapply(1:10, FUN=function(i) {
    print(RF_smote_sampling$bestTune$mtry)
    RF_smote_sampling <-
        balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                           y_train_factor = y_all_data,
                                           model_method = "rf",
                                           fitControl = fitControl,
                                           sampling_method = "smote",
                                           tuneGrid =
                                               data.frame(mtry=RF_smote_sampling$bestTune$mtry),
                                           verboseIter = TRUE)
    var.imp.rf <- varImp(RF_smote_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance),
                                   var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)),
                                   decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
})

rf_smote_sampling_selected_genes <- sort(table(unlist(lapply(
    RF_smote_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]
RF_smote_selected_genes_eczI_lymnI <- gsub("X","",names(rf_smote_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------
# nI I
d_train <- rbind(v_ecz_nI, v_lym_I)
y_train_factor <- as.factor(d_train$class)
d_train$class <- NULL


RF_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                        y_train_factor = y_train_factor,
                                                        model_method = "rf",
                                                        fitControl = fitControl,
                                                        sampling_method = "smote",verboseIter = TRUE)


RF_smote_sampling_selected_genes <- lapply(1:10, FUN=function(i) {
    print(RF_smote_sampling$bestTune$mtry)
    RF_smote_sampling <-
        balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                           y_train_factor = y_all_data,
                                           model_method = "rf",
                                           fitControl = fitControl,
                                           sampling_method = "smote",
                                           tuneGrid =
                                               data.frame(mtry=RF_smote_sampling$bestTune$mtry),
                                           verboseIter = TRUE)
    var.imp.rf <- varImp(RF_smote_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance),
                                   var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)),
                                   decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
})

rf_smote_sampling_selected_genes <- sort(table(unlist(lapply(
    RF_smote_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]
RF_smote_selected_genes_ecznI_lymI <- gsub("X","",names(rf_smote_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------
# nI nI
d_train <- rbind(v_ecz_I, v_lym_I)
y_train_factor <- as.factor(d_train$class)
d_train$class <- NULL


RF_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                        y_train_factor = y_train_factor,
                                                        model_method = "rf",
                                                        fitControl = fitControl,
                                                        sampling_method = "smote",verboseIter = TRUE)


RF_smote_sampling_selected_genes <- lapply(1:10, FUN=function(i) {
    print(RF_smote_sampling$bestTune$mtry)
    RF_smote_sampling <-
        balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                           y_train_factor = y_all_data,
                                           model_method = "rf",
                                           fitControl = fitControl,
                                           sampling_method = "smote",
                                           tuneGrid =
                                               data.frame(mtry=RF_smote_sampling$bestTune$mtry),
                                           verboseIter = TRUE)
    var.imp.rf <- varImp(RF_smote_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance),
                                   var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)),
                                   decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
})

rf_smote_sampling_selected_genes <- sort(table(unlist(lapply(
    RF_smote_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]
RF_smote_selected_genes_I_I <- gsub("X","",names(rf_smote_sampling_selected_genes)[-1])



i <-intersect(intersect(RF_smote_selected_genes_eczI_lymnI,RF_smote_selected_genes_ecznI_lymI),
    intersect(RF_smote_selected_genes_I_I,RF_smote_selected_genes_nI_nI))


#--------------------------------------------------------------------
MDS_for_selected_genes(dds.ecz = dds.ecz,dds.lym = dds.lym,
                       sel_genes = list(i,
                                        RF_smote_selected_genes_eczI_lymnI,
                                        RF_smote_selected_genes_ecznI_lymI,
                                        RF_smote_selected_genes_I_I,
                                        RF_smote_selected_genes_nI_nI),
                       titles = c("intersect",
                                  "RF_smote_selected_genes_eczI_lymnI",
                                  "RF_smote_selected_genes_ecznI_lymI",
                                  "RF_smote_selected_genes_I_I",
                                  "RF_smote_selected_genes_nI_nI"),
                       include_only_disease = FALSE,measurement = "rpkm",ebg = ebg)


compare_selected_features_barplot(features = i,dds.lym = dds.lym,
                                  lym.norm = lym.norm,ecz.norm = ecz.norm,
                                  dds.ecz = dds.ecz)
