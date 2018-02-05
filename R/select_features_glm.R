#-------------------------------------------------------------------------------
# feature selection with lasso
# create a model.matrix of train  and test features, subseting just removing the intercept
# x_train <-  model.matrix(class~., d_with_class_train) [ , -1]
# x_test <- model.matrix(class~., d_with_class_test) [ , -1]
#
# # lasso model with glmnet
# lasso_mod_cv = cv.glmnet(x_train, y_train, alpha =1)
# bestlam = lasso_mod_cv$lambda.min
# predict(lasso_mod_cv,newdata = x_test,type = "coefficients", s = bestlam )
#
# p <- ifelse(predict(lasso_mod_cv ,type ="response", s = bestlam , newx =x_test) <.5, 0,1)
# table(y_test,p)

#-------------------------------------------------------------------------------
# lasso with caret + glmnet
# Tuning parameters
# 1. alpha
# 2. lambda
# train the model
# tuneGrid=expand.grid(.alpha=1,.lambda=seq(0, 100, by = 0.1))  # tuneGrid = tuneGrid,

# train one model to find the best tuning parameters
glmFit.up <- train(x = d_train, y = y_train_factor,
                   method = "glmnet",
                   trControl = fitControl,
                   metric = "ROC")

# train the model several times using only the best parameters
glm.feats <- lapply(1:20, FUN=function(i) {
    cat(i)
    glmFit.up <- train(x = d_train, y = y_train_factor,
                       method = "glmnet",
                       trControl = fitControl,
                       tuneGrid = data.frame(alpha = glmFit.up$bestTune$alpha,
                                            lambda = glmFit.up$bestTune$lambda),
                       metric = "ROC")
    var.imp.glm <- varImp(glmFit.up) # glmnet

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




glm <- sort(table(unlist(lapply(glm.feats, FUN=function(x)x$genes))), decreasing=T)[1:15]
glm <- cbind(names(glm), glm)
#

