glmnet_original<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                y_train_factor = y_train_factor,
                                                model_method = "glmnet",
                                                fitControl = fitControl,
                                                sampling_method = "original",
                                                family = "binomial")

glmnet_down_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                     y_train_factor = y_train_factor,
                                                     model_method = "glmnet",
                                                     fitControl = fitControl,
                                                     sampling_method = "down",
                                                     family = "binomial")

glmnet_up_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                   y_train_factor = y_train_factor,
                                                   model_method = "glmnet",
                                                   fitControl = fitControl,
                                                   sampling_method = "up",
                                                   family = "binomial")

glmnet_smote_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                      y_train_factor = y_train_factor,
                                                      model_method = "glmnet",
                                                      fitControl = fitControl,
                                                      sampling_method = "smote",
                                                      family = "binomial")

glmnet_weighted_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                         y_train_factor = y_train_factor,
                                                         model_method = "glmnet",
                                                         fitControl = fitControl,
                                                         sampling_method = "weights",
                                                         family = "binomial")



model_list = list("glmnet_original" = glmnet_original,
                  "glmnet_up_sampling" = glmnet_up_sampling,
                  "glmnet_down_sampling" = glmnet_down_sampling,
                  "glmnet_weighted_sampling" = glmnet_weighted_sampling,
                  "glmnet_smote_sampling" = glmnet_smote_sampling)
plot_multipe_rocs(model_list,
                  use_test_data = FALSE,
                  y_test_factor = y_test_factor,
                  d_with_class_test = d_with_class_test,
                  variable_of_interest = "ecz")
compare_sampling_approaches_on_test_data(model_list,d_with_class_test)
