# loading libraries
library(dplyr) # for data manipulation
library(caret) # for model-building
library(DMwR) # for smote implementation
library(purrr) # for functional programming (map)
library(pROC) # for AUC calculations
#-------------------------------------------------------------------------------
# configure caret
# parameter tuning
library(doMC)
registerDoMC(cores = 3)
# Create a stratified sample for repeated cv
#cv_10_folds <- createMultiFolds(y_train_factor, k=10, times =2)
#index = cv_10_folds
fitControl <- trainControl(method = "cv",
                           classProbs = T,
                           savePredictions = T,
                           summaryFunction = twoClassSummary,
                           verboseIter = TRUE, allowParallel = TRUE)

