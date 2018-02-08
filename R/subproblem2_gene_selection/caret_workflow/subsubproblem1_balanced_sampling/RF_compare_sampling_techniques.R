RF_original<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                     y_train_factor = y_train_factor,
                                                     model_method = "rf",
                                                     fitControl = fitControl,
                                                     sampling_method = "original")

RF_down_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                     y_train_factor = y_train_factor,
                                                     model_method = "rf",
                                                     fitControl = fitControl,
                                                     sampling_method = "down")

RF_up_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                   y_train_factor = y_train_factor,
                                                   model_method = "rf",
                                                   fitControl = fitControl,
                                                   sampling_method = "up")

RF_smote_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                   y_train_factor = y_train_factor,
                                                   model_method = "rf",
                                                   fitControl = fitControl,
                                                   sampling_method = "smote")

RF_weighted_sampling<-balanced_sampling_methods_with_ROC(d_train = d_train,
                                                   y_train_factor = y_train_factor,
                                                   model_method = "rf",
                                                   fitControl = fitControl,
                                                   sampling_method = "weights")



model_list = list("RF_original" = RF_original,
                  "RF_up_sampling" = RF_up_sampling,
                  "RF_down_sampling" = RF_down_sampling,
                  "RF_weighted_sampling" = RF_weighted_sampling,
                  "RF_smote_sampling" = RF_smote_sampling)
plot_multipe_rocs(model_list,
                  use_test_data = FALSE,
                  y_test_factor = y_test_factor,
                  d_with_class_test = d_with_class_test,
                  variable_of_interest = "ecz")
compare_sampling_approaches_on_test_data(model_list,d_with_class_test)
