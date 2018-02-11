glm_up_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                      y_train_factor = y_train_factor,
                                                      model_method = "glmnet",
                                                      fitControl = fitControl,
                                                      sampling_method = "up",
                                                      family = "binomial")

glm_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                      y_train_factor = y_train_factor,
                                                      model_method = "glmnet",
                                                      fitControl = fitControl,
                                                      sampling_method = "smote",
                                                      family = "binomial")

RF_up_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                     y_train_factor = y_train_factor,
                                                     model_method = "rf",
                                                     fitControl = fitControl,
                                                     sampling_method = "up",verboseIter = TRUE)


RF_smote_sampling <- balanced_sampling_methods_with_ROC(d_train = d_train,
                                                     y_train_factor = y_train_factor,
                                                     model_method = "rf",
                                                     fitControl = fitControl,
                                                     sampling_method = "smote",verboseIter = TRUE)



#-------------------------------------------------------------------------------
model_list = list(RF_up_sampling,
                  RF_smote_sampling,
                  glm_up_sampling,
                  glm_smote_sampling)

titles = c("RF_up_sampling",
           "RF_smote_sampling",
           "glm_up_sampling",
           "glm_smote_sampling")

compare_sampling_approaches_on_test_data(model_list, titles,
                                         d_with_class_test = d_with_class_test)
#-------------------------------------------------------------------------------
# compare ROC
plot_multipe_rocs(model_list, titles)

