# random forest with ebm

# [question] what exactly is the tuning parameter here ?
# step 1. build a model
rfFit.up <- caret::train(x = d_train, y = y_train_factor,
                         method="rf",
                         trControl =fitControl,
                         metric = "ROC")


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



