# Build a standard classifier using a gradient boosted machine
# unbalanced sampling - didn't work BTW
orig_fit <- train(y = y_train_factor,
                  x = d_with_class_train,
                  method = "gbm",
                  verbose = FALSE,
                  metric = "ROC",
                  trControl = fitControl)




#-------------
ifelse(predict(weighted_fit, d_test, type = "prob")[,"ecz"]>.5,0,1) %>%
    confusionMatrix(y_test)
r<-test_roc(weighted_fit, d_with_class_test)
