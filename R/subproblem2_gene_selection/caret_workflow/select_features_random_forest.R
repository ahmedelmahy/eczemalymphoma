# random forest with ebm

# [question] what exactly is the tuning parameter here ?
# step 1. build a model
rfFit <- caret::train(x = d_train, y = y_train_factor,
                         method="rf",
                         trControl =fitControl,
                         metric = "ROC")

plot(varImp(train_cv_1), main = "Random Forest - Variable Importance plot")

l <- list(NA)
a<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                   y_train_factor = y_train_factor,
                                   d_test = d_test,
                                   y_test_factor = y_test_factor,
                                   model_method = "rf",
                                   fitControl = fitControl,
                                   y_factor_to_predict = y_factor_to_predict,
                                   sampling_method = "smote",
                                   d_with_class_test = d_with_class_test)



b<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                         y_train_factor = y_train_factor,
                                         d_test = d_test,
                                         y_test_factor = y_test_factor,
                                         model_method = "rf",
                                         fitControl = fitControl,
                                         y_factor_to_predict = y_factor_to_predict,
                                         sampling_method = "up",
                                         d_with_class_test = d_with_class_test)

plot_multipe_rocs(list(a,b))
length(l)

l[1] <- 4
l
rf.feats <- lapply(1:20, FUN=function(i) {
    cat(i)
    rfFit.up <- caret::train(x = d_train, y = y_train_factor,
                             method="rf",
                             trControl = fitControl,
                             tuneGrid = data.frame(mtry=rfFit.up$bestTune$mtry),
                             metric = "ROC")
    var.imp.rf <- varImp(rfFit.up) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance), var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)), decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
})

rf <- sort(table(unlist(lapply(rf.feats, FUN=function(x)x$genes))), decreasing=T)[1:15]
rf <- cbind(names(rf), rf)



