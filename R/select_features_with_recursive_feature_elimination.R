rf.vs <- varSelRF(data.mat, class.labels) # using oob error for recurssive feature elimination
oobrf <- rf.vs$selected.vars
### parameter tuning with selected features
sel.genes <- c(9760, 8631)
sel.data <- new.data[, c(paste0("X", sel.genes), "class")]




# ! ?
sel.data <- data.frame(sel.data)
sel.data$class <- class.labels
