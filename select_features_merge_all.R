sel.features <- merge(glm, rf, all=T)
#sel.features <- as.matrix(merge(sel.features, lm, all=T))

sel.features[which(is.na(sel.features))] <- 0
sel.features <- data.frame(sel.features)


# pdf("sel.features_up.pdf")
# pheatmap(data.matrix(sel.features[, 2:4]), cluster_rows = T, cluster_cols = F, labels_row = sel.features$V1, legend = F, cellwidth = 20, scale = "column")
# pheatmap(data.matrix(sel.features[, 2:4]), cluster_rows = T, cluster_cols = F, labels_row = paste0("gene", 1:nrow(sel.features)), legend = F, cellwidth = 35, scale = "column")
# dev.off()

#-------------------------------------------------------------------------------
# implement stepwise variable selection on the sel.features
f <- formula(paste("class ~ ",paste(sel.features[[1]],collapse = " + ")))
# full model with all variables
fullmod = glm(f, data = d_with_class_train,family=binomial)
# forward and backward variable selection
backwards = step(fullmod, direction = "both") # selects only three variables
## same result from
## library(MASS)
## step = stepAIC(fullmod, direction = "both", trace = 0)

# full glm on lasso variables only

f <- formula(paste("class ~ ",paste(rownames(glm),collapse = " + ")))
fullmod_glm = glm(f, data = d_with_class_train,family=binomial)

#-------------- perfom lasso on the orginal data to compare > well it out perfoms them!
# but if we only did one lasso
x_train <-  model.matrix(class~., d_with_class_train) [ , -1]
x_test <- model.matrix(class~., d_with_class_test) [ , -1]

lasso_mod_cv = cv.glmnet(x_train, y_train, alpha =1)
bestlam = lasso_mod_cv$lambda.min
#predict(lasso_mod_cv,newdata = x_test,type = "coefficients", s = bestlam )



# the three variables perform better on the test but all are bad
p_full_glm <- ifelse(predict(fullmod_glm, newdata = d_test, type="response") < .5,0,1)
p_full <- ifelse(predict(fullmod, newdata = d_test, type="response") < .5,0,1)
p_step <- ifelse(predict(backwards, newdata = d_test, type="response") < .5,0,1)
p_lasso <- ifelse(predict(lasso_mod_cv ,type ="response", s = bestlam , newx =x_test) <.5, 0,1)

table(y_test, p_step)
table(y_test, p_full)
table(y_test, p_full_glm)
table(y_test,p_lasso)

