# svm
svmFit.up <- train(class ~ ., data = new.data,
                   method = "svmRadial",
                   trControl = fitControl,
                   preProc = c("center", "scale"),
                   tuneLength = 8,
                   metric = "ROC")


# linear svm
lsvmFit.up <- train(class ~ ., data = new.data,
                    method = "svmLinear",
                    trControl = fitControl,
                    preProc = c("center", "scale"),
                    tuneLength = 8,
                    metric = "ROC")
