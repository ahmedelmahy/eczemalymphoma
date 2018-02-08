# selecting important variables using random forest
rf.feats <- lapply(1:2, FUN=function(i) {
    if (i == 1){
        RF_up_sampling <-
            balanced_sampling_methods_with_ROC(d_train = d_train,
                                               y_train_factor = y_train_factor,
                                               model_method = "rf",
                                               fitControl = fitControl,
                                               sampling_method = "up")
        i <- 2
    }
    if (i > 1) {
    RF_up_sampling <-
        balanced_sampling_methods_with_ROC(d_train = d_train,
                                                       y_train_factor = y_train_factor,
                                                       model_method = "rf",
                                                       fitControl = fitControl,
                                                       sampling_method = "up",
                                                       tuneGrid =
                                               data.frame(mtry=RF_up_sampling$bestTune$mtry))
    var.imp.rf <- varImp(RF_up_sampling) # random forest
    var.imp.rf <- data.frame(cbind(row.names(var.imp.rf$importance), var.imp.rf$importance[, 1]))
    names(var.imp.rf) <- c("genes", "imp")
    var.imp.rf <- var.imp.rf[order(as.numeric(as.matrix(var.imp.rf$imp)), decreasing = T), ]
    #var.imp.rf$genes <- sub("X", "", var.imp.rf$genes)
    return(var.imp.rf[1:20, ])
    }
})

rf <- sort(table(unlist(lapply(rf.feats, FUN=function(x)x$genes))), decreasing=T)

rf <- data.frame(rf)



