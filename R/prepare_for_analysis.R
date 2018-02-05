set.seed(85)
# load additional functions
source("R/classifier_funcs.R")
source("R/load.R")
source("R/read_lymphoma_data.R")
source("R/read_eczema_data.R")
source("R/merge_eczema_lymphoma_data.R")

#register multicore parallel
library(TSA)
#library(doMC)
#registerDoMC(2)

#prepare the data
# remove X from gene names
#colnames(new.data)<- sub("X", "", colnames(new.data))  # I guess didn't work
# because variables as only numbers doesn't fit in formulas
#-------------------------------------------------------------------------------
# step 0 tidy the data
# two datasets with only numbers

# binarize the class variable
new.data$class <- ifelse(new.data$class == "lym",1,0)  # models will predict lymphoma

#-------------------------------------------------------------------------------
# split the data into train and test to compare different models
train= sample(dim(new.data)[1],30)
# for models that require class to be part of the dataset
d_with_class_train <- new.data[train,]
d_with_class_test <- new.data[-train,]

# for models that require class to be a separate variable
# and for models that require the class to be factor
d_train <- d_with_class_train
y_train <- d_train$class
y_train_factor <- as.factor(ifelse(y_train == 0, "ecz","lym"))
d_train$class <- NULL

d_test <- d_with_class_test
y_test <- d_test$class
y_test_factor <- as.factor(ifelse(y_test == 0, "ecz","lym"))

d_test$class <- NULL



#-------------------------------------------------------------------------------
# configure caret
# parameter tuning
# modify the methods used by Caret for CV
# fitControl(list)
fitControl <- trainControl(method = "cv", sampling = "up",
                           classProbs = T, summaryFunction = twoClassSummary,
                           verboseIter = TRUE)

#fitControl <- trainControl(method = "none",classProbs = TRUE,summaryFunction = twoClassSummary,
 #                          verboseIter = TRUE)

#fitControl <- trainControl(method = "LOOCV", sampling = "up", classProbs = T, summaryFunction = twoClassSummary)
