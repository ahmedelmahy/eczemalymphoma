
test_roc <- function(model, data) {
    roc(data$class,
        predict(model, data, type = "prob")[, "ecz"])
}


balanced_sampling_methods_with_ROC <- function(d_train,
                                               y_train_factor,
                                               d_test,
                                               y_test_factor,
                                               model_method,
                                               fitControl,
                                               y_factor_to_predict,
                                               sampling_method,
                                               d_with_class_test){
    switch(sampling_method,
           original = {
               message("doing the original model")
               fit <- caret::train(x = d_train, y = y_train_factor,
                                   method=model_method,
                                   trControl =fitControl,
                                   metric = "ROC")

               message("#fitControl now uses the same seed to ensure same cross-validation splits")
               fitControl$seeds <<- fit$control$seeds
           } ,
           weights = {
               message("#using weights for sampling")
               model_weights <- ifelse(y_train_factor ==
                                           names(table(y_train_factor)[1]),
                                       (1/table(y_train_factor)[1]) * 0.5,
                                       (1/table(y_train_factor)[2]) * 0.5)

               fit <- caret::train(x = d_train, y = y_train_factor,
                                   method=model_method,
                                   weights = model_weights,
                                   trControl =fitControl,
                                   metric = "ROC")
           },
           up = {
               message("# using up sampling")
               fitControl$sampling <- "up"
               fit <- train(y = y_train_factor,
                            x = d_with_class_train,
                            method =model_method,
                            verbose = FALSE,
                            metric = "ROC",
                            trControl = fitControl)
           },
           down = {
               message("# using down sampling")
               fitControl$sampling <- "down"
               fit <- train(y = y_train_factor,
                            x = d_with_class_train,
                            method = model_method,
                            verbose = FALSE,
                            metric = "ROC",
                            trControl = fitControl)
           },
           smote = {
               message("# using smote sampling")
               fitControl$sampling <- "smote"
               fit <- train(y = y_train_factor,
                            x = d_with_class_train,
                            method = model_method,
                            verbose = FALSE,
                            metric = "ROC",
                            trControl = fitControl)
           })
    fit_roc <- test_roc(model = fit, data = d_with_class_test)
    message(auc(fit_roc))
    return(fit_roc)
}

plot_multipe_rocs <- function(model_list_roc){
    results_list_roc <- list(NA)
    num_mod <- 1
    for(the_roc in model_list_roc){

        results_list_roc[[num_mod]] <-
            data_frame(tpr = the_roc$sensitivities,
                       fpr = 1 - the_roc$specificities,
                       model = names(model_list)[num_mod])
        num_mod <- num_mod + 1
    }
    results_df_roc <- bind_rows(results_list_roc)
    # Plot ROC curve for all 5 models
    ggplot(aes(x = fpr,  y = tpr, group = model), data = results_df_roc) +
        geom_line(aes(color = model), size = 1) +
        geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
        theme_bw(base_size = 18)
}




