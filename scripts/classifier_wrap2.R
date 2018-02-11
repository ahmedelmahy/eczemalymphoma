# random forest
rfFit.up <- caret::train(class ~ ., data = sel.data,
                         method="rf",
                         trControl = fitControl,
                         metric = "ROC")
# svm
svmFit.up <- train(class ~ ., data = sel.data,
                   method = "svmRadial",
                   trControl = fitControl,
                   preProc = c("center", "scale"),
                   tuneLength = 8,
                   metric = "ROC")
# linear svm
lsvmFit.up <- train(class ~ ., data = sel.data,
                    method = "svmLinear",
                    trControl = fitControl,
                    preProc = c("center", "scale"),
                    tuneLength = 8,
                    metric = "ROC")
# glm
glmFit.up <- train(class ~ ., data = sel.data,
                   method = "glmnet",
                   trControl = fitControl)

### final model with tuned parameters and selected features
test.svms <- lapply(1:100, FUN=function(i){
    test.svm <- train(class ~ ., data = sel.data,
                      method = "svmRadial",
                      trControl = fitControl,
                      tuneGrid = data.frame(sigma = svmFit.up$bestTune$sigma, C=svmFit.up$bestTune$C),
                      metric = "ROC")
    return(test.svm$results)
})

test.lsvms <- lapply(1:100, FUN=function(i){
    test.lsvm <- train(class ~ ., data = sel.data,
                       method = "svmLinear",
                       trControl = fitControl,
                       tuneGrid = data.frame(C=lsvmFit.up$bestTune$C),
                       metric = "ROC")
    return(test.lsvm$results)
})
test.glms <- lapply(1:100, FUN=function(i){
    test.glm.up <- train(class ~ ., data = sel.data,
                         method = "glmnet",
                         trControl = fitControl,
                         tuneGrid = data.frame(alpha = glmFit.up$bestTune$alpha, lambda = glmFit.up$bestTune$lambda),
                         metric = "ROC")
    return(test.glm.up$results)
})

test.rfs <- lapply(1:100, FUN=function(i){
    test.rf.up <- caret::train(class ~ ., data = sel.data,
                               method="rf",
                               trControl = fitControl,
                               tuneGrid = data.frame(mtry=rfFit.up$bestTune$mtry),
                               metric = "ROC")
    return(test.rf.up$results)
})

# roc plots

glm.rocs <- data.frame(do.call(rbind, test.glms))
glm.rocs$model <- "glm"
rf.rocs <- data.frame(do.call(rbind, test.rfs))
rf.rocs$model <- "rf"
svm.rocs <- data.frame(do.call(rbind, test.svms))
svm.rocs$model <- "svm"
lsvm.rocs <- data.frame(do.call(rbind, test.lsvms))
lsvm.rocs$model <- "lsvm"

all.rocs <- rbind.fill(list(glm.rocs, rf.rocs, svm.rocs, lsvm.rocs))
pdf("rocs_up2.pdf")
boxplot(split(all.rocs$ROC, all.rocs$model), ylim=c(0,1))

plot.mat <- melt(all.rocs, id=c("alpha", "lambda", "mtry", "sigma", "C", "model"))
# model comparison
levels(plot.mat$variable)[levels(plot.mat$variable)=="Sens"] <- "Sensitivity"
levels(plot.mat$variable)[levels(plot.mat$variable)=="Spec"] <- "Specificity"

p <- ggplot(plot.mat, aes(y=plot.mat$value, x=plot.mat$model)) +
    facet_grid(variable ~ .) + geom_boxplot() + ylim(0,1) +
    labs(title="") + xlab("") + ylab("") +
    theme_bw() + theme(legend.title=element_blank())+
    theme(legend.position="none")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
print(p)
dev.off()
