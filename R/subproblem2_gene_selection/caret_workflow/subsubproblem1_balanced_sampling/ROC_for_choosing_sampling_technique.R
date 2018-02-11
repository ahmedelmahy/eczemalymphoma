balanced_sampling_methods_with_ROC <- function(d_train,
                                               y_train_factor,
                                               model_method,
                                               fitControl,
                                               sampling_method, ...){
    switch(sampling_method,
           original = {
               message("doing the original model")
               fit <- caret::train(x = d_train, y = y_train_factor,
                                   method=model_method,
                                   trControl =fitControl,
                                   metric = "ROC",...)

               #message("#fitControl now uses the same seed to ensure same cross-validation splits")
               #fitControl$seeds <<- fit$control$seeds
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
                                   metric = "ROC", ...)
           },
           up = {
               message("# using up sampling")
               fitControl$sampling <- "up"
               fit <- train(x = d_train,
                            y = y_train_factor,
                            method =model_method,
                            metric = "ROC",
                            trControl = fitControl,...)

           },
           down = {
               message("# using down sampling")
               fitControl$sampling <- "down"
               fit <- train(x = d_train,
                            y = y_train_factor,
                            method = model_method,
                            metric = "ROC",
                            trControl = fitControl,...)
           },
           smote = {
               message("# using smote sampling")
               fitControl$sampling <- "smote"
               fit <- train(x = d_train,
                            y = y_train_factor,
                            method = model_method,
                            metric = "ROC",
                            trControl = fitControl,...)
           })
    return(fit)
}

plot_multipe_rocs <- function(model_list,titles, use_test_data = FALSE,
                              y_test_factor = NA, d_with_class_test = NA,
                              variable_of_interest = "ecz"){
    if (length(model_list) != length(titles)) {stop("titles are not the same
                                                       length as selected_genes_list")}
    if (use_test_data){
        results_list_roc <- list(NA)
        num_mod <- 1
        # fit = model_list[1]
        for(fit in model_list){
            message("ROC for test data")
            fit_probs <- predict(fit, d_with_class_test,
                                 type = "prob")[, variable_of_interest]
            #message(auc(fit_roc))
            # roc for the test dataset
            myRoc_on_test_data <- roc(predictor=fit_probs,
                                      response=y_test_factor,
                                      levels=rev(levels(y_test_factor)))

            message(auc(myRoc_on_test_data))
            results_list_roc[[num_mod]] <- data_frame(
                True_positive_rate = myRoc_on_test_data$sensitivities,
                False_positive_rate = 1 - myRoc_on_test_data$specificities,
                model = titles[num_mod])
            num_mod <- num_mod + 1 }
    } else {
        message("ROC for training data")
        results_list_roc <- list(NA)
        num_mod <- 1
        for(model in model_list){
            myRoc_on_train_data <- roc(predictor = model$pred$lym,
                                       response = model$pred$obs,
                                       positive = variable_of_interest)
            message(auc(myRoc_on_train_data))
            results_list_roc[[num_mod]] <- data_frame(
                True_positive_rate = myRoc_on_train_data$sensitivities,
                False_positive_rate = 1 - myRoc_on_train_data$specificities,
                model = titles[num_mod])
            num_mod <- num_mod + 1
        }
    }

    results_df_roc <- bind_rows(results_list_roc)
    # Plot ROC curve for all 5 models
    ggplot(aes(x = False_positive_rate,  y = True_positive_rate, group = model), data = results_df_roc) +
        geom_line(aes(color = model), size = 1) +
        geom_abline(intercept = 0, slope = 1, color = "gray", size = 1) +
        theme_bw(base_size = 18)
}

# compare perfomance of different sampling approaches on test set
compare_sampling_approaches_on_test_data <- function(model_list,
                                                     titles,
                                                     d_with_class_test){
    if (length(model_list) != length(titles)) {stop("titles are not the same
                                                       length as selected_genes_list")}
    num_mod <- 1
    df <- tibble("sampling_approach" = 0, "class" = "0",
                 "probability_of_success" = 0)
    #fit <- RF_up_sampling
for(fit in model_list){
    fit_probs <- predict(fit, d_with_class_test,type = "prob")
    l <- dim(fit_probs)[1]
    for (i in 1:l) {
        probability_of_success <- fit_probs[i,y_test_factor[i]]
        df[dim(df)[1]+1,] <- data_frame(
            "sampling_approach" = titles[num_mod],
            "class" = as.character(y_test_factor[i]),
            "probability_of_success" = probability_of_success
        )
    }
    num_mod <- num_mod + 1
}
    df<- df[-1,]
    df_grouped <-df %>%
        group_by(sampling_approach,class)%>%
        summarise(mean(probability_of_success)) %>%
        set_colnames(c("sampling_approach","class",
                 "mean_of_probabilities_of_choosing_the_correct_class"))
    p <-ggplot(df_grouped,aes(sampling_approach,
                          mean_of_probabilities_of_choosing_the_correct_class,
                          fill = class))+
        geom_bar(stat="identity", position=position_dodge())+
        theme(axis.text.x = element_text(angle = 90))

    print(p )



    }


