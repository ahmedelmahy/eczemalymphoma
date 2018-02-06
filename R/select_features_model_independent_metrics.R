# [question] didn't get it
# Model Independent Metrics
# filterVarImp calculates the area under the ROC curve for every variable
# "independently" for each class.
# this is the defaut way to get scores from models where no built-in
# importance score is implemented (or exists) such as svm
var.imp <- filterVarImp(d_train, y_train_factor) # linear model

var.imp <- var.imp[order(var.imp$ecz, decreasing = T), ]
#rownames(var.imp) <- sub("X", "", rownames(var.imp))
lm <- cbind(rownames(var.imp[1:10, ]), round(var.imp[1:10, 1]*100))
