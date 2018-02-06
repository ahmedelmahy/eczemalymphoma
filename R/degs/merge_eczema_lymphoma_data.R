# genes to be used for classifier
# [question] genes that are in both datasets (the differentially expressed genes that maybe in either dataset or gene.mat)
gene.mat <- merge(data.frame(lym.sel), data.frame(ecz.sel), by="entrez", all=T)
if(TRUE) cat("\n Dimentions of lymphoma genes ",dim(data.frame(lym.sel))," \n\t eczema genes",
    dim(data.frame(ecz.sel)),"\n\t now all genes included in gene mat are",dim(gene.mat))
# data matrix from the two classes
lym.mat <- data.frame(lym.fcs[which(rownames(lym.fcs)%in%gene.mat$entrez), ])
ecz.mat <- data.frame(ecz.fcs[which(rownames(ecz.fcs)%in%gene.mat$entrez), ])
data.mat <- cbind(lym.mat, ecz.mat[rownames(lym.mat), ])


### selected features matrix
# remove NA and NANs ; already removed before so 0 is removed
# and transpose
data.mat <- t(data.mat[which(apply(data.mat, 1, function(x) length(which(is.infinite(x) | is.na(x))))==0), ])
# create labels for classification either lymphoma or eczema
class.labels <- as.factor(c(rep("lym", ncol(lym.mat)), rep("ecz", ncol(ecz.mat))))

# change data.mat to a dataframe and add the classification labels
new.data <- data.frame(data.mat)
new.data$class <- class.labels

# now we have a dataset ready for analysis
#---------------------------------------------------------------------------------
