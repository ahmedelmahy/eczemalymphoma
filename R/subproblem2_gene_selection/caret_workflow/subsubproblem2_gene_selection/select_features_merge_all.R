sel.features <- merge(glm, rf, all=T)
#sel.features <- as.matrix(merge(sel.features, lm, all=T))

sel.features[which(is.na(sel.features))] <- 0
sel.features <- data.frame(sel.features)



sel.data <- new.data[, c(paste0("X", sel.genes), "class")]
# ! ?
sel.data <- data.frame(sel.data)
sel.data$class <- class.labels

# pdf("sel.features_up.pdf")
# pheatmap(data.matrix(sel.features[, 2:4]), cluster_rows = T, cluster_cols = F, labels_row = sel.features$V1, legend = F, cellwidth = 20, scale = "column")
# pheatmap(data.matrix(sel.features[, 2:4]), cluster_rows = T, cluster_cols = F, labels_row = paste0("gene", 1:nrow(sel.features)), legend = F, cellwidth = 35, scale = "column")
# dev.off()


# full glm on lasso variables only

f <- formula(paste("class ~ ",paste(rownames(glm),collapse = " + ")))
fullmod_glm = glm(f, data = d_with_class_train,family=binomial)

#-------------- perfom lasso on the orginal data to compare > well it out perfoms them!
# but if we only did one lasso



# the three variables perform better on the test but all are bad
p_full_glm <- ifelse(predict(fullmod_glm, newdata = d_test, type="response") < .5,0,1)
p_full <- ifelse(predict(fullmod, newdata = d_test, type="response") < .5,0,1)


table(y_test, p_full)
table(y_test, p_full_glm)


