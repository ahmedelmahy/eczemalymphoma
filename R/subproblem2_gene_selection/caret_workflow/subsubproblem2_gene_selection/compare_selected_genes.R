RF_smote_sampling_selected_genes <- lapply(1:2, FUN=function(i) {
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
RF_smote_selected_genes <- gsub("X","",names(rf_smote_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------
RF_up_sampling_selected_genes <- lapply(1:2, FUN=function(i) {
    print(RF_up_sampling$bestTune$mtry)
    RF_up_sampling <-
        balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                           y_train_factor = y_all_data,
                                           model_method = "rf",
                                           fitControl = fitControl,
                                           sampling_method = "up",
                                           tuneGrid =
                                               data.frame(mtry=RF_up_sampling$bestTune$mtry),
                                           verboseIter = TRUE)
    var.imp.rf <- varImp(RF_up_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance),
                                   var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)),
                                   decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
    })

rf_up_sampling_selected_genes <- sort(table(unlist(lapply(
    RF_up_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]
RF_up_selected_genes <- gsub("X","",names(rf_up_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------
glm_up_sampling_selected_genes <- lapply(1:2, FUN=function(i) {
    cat(i)
    glm_up_sampling <- balanced_sampling_methods_with_ROC(d_train = all_data_genes,
                                                          y_train_factor = y_all_data,
                                                          model_method = "glmnet",
                                                          fitControl = fitControl,
                                                          sampling_method = "up",
                                                          family = "binomial",
                                                          tuneGrid = data.frame(
                                                              alpha = glm_up_sampling$bestTune$alpha,
                                                              lambda = glm_up_sampling$bestTune$lambda))
    var.imp.glm <- varImp(glm_up_sampling) # glmnet

    # find top n important variables
    var.imp.glm <- data.frame(genes = row.names(var.imp.glm$importance),
                              imp = var.imp.glm$importance[, 1])
    var.imp.glm <- var.imp.glm[order(var.imp.glm$imp, decreasing = T), ]
    #[question] well, top 20 , but only first 3 has values
    #  so remove genes with 0 coefficient, or instead of that we remove NAs from
    #  the returned df
    #var.imp.glm <- var.imp.glm[var.imp.glm$imp != 0,]

    return(na.omit(var.imp.glm[1:20, ]))
})
glm_up_sampling_selected_genes <- sort(table(unlist(lapply(
    glm_up_sampling_selected_genes,FUN=function(x)x$genes))), decreasing=T)[1:30]

glm_up_selected_genes <- gsub("X","",names(glm_up_sampling_selected_genes)[-1])

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# load more selected variables
#source("R/subproblem1_degs/method1_DESeq/appraoch2_disease1normal2disease2normal1/multiple_comparision_then_intersection.R")
#-------------------------------------------------------------------------------
# compare them

i <- intersect(eIlI,glm_up_selected_genes)
i <- intersect(eIlI, RF_up_selected_genes)
i <- intersect(glm_up_selected_genes, RF_up_selected_genes)
i <- intersect(i, eIlI)
i <- intersect(setdiff(glm_up_selected_genes,enIlnI),
               setdiff(RF_up_selected_genes,enIlnI))

i <- intersect(eIlI,enIlnI)

glm_intersect_enIlI_u_eIlnI <- intersect(glm_up_selected_genes,union(enIlI,eIlnI))
glm_i_enIlI_u_eIlnI_i_eIlI <- intersect(intersect(glm_up_selected_genes,
                                            union(enIlI,eIlnI)),eIlnI)

RF_i_enIlI_u_eIlnI <- intersect(RF_up_selected_genes,union(enIlI,eIlnI))
RF_d_enIlI_i_eIlnI_d_enIlnI <- setdiff(intersect(RF_up_selected_genes,
                                               union(enIlI,eIlnI)),enIlnI)

glm_d_enIlI_u_eIlnI <- setdiff(glm_up_selected_genes,union(enIlI,eIlnI))
most_genes <- sample(union(union(eIlI,enIlnI),union(enIlI,eIlnI)),100)
MDS_for_selected_genes(dds.ecz = dds.ecz,dds.lym = dds.lym,
                       sel_genes = list(most_genes,
                                        RF_up_selected_genes,
                                        glm_up_selected_genes,
                                        eIlI,
                                        RF_i_enIlI_u_eIlnI,
                                        RF_d_enIlI_i_eIlnI_d_enIlnI,
                                        glm_i_enIlI_u_eIlnI_i_eIlI,
                                        RF_smote_selected_genes),
                       titles = c("random 100 genes",
                                  "RF_up_selected_genes",
                                  "glm_up_selected_genes",
                                  "eIlI",
                                  "RF_i_enIlI_u_eIlnI",
                                  "RF_d_enIlI_i_eIlnI_d_enIlnI",
                                  "glm_i_enIlI_u_eIlnI_i_eIlI",
                                  "RF_smote_selected_genes"),
                       include_only_disease = TRUE,measurement = "rpkm",ebg = ebg)


i <-intersect(RF_smote_selected_genes,RF_up_selected_genes)
compare_selected_features_barplot(features = i,dds.lym = dds.lym,
                                  lym.norm = lym.norm,ecz.norm = ecz.norm,
                                  dds.ecz = dds.ecz)

