# loading required packages
# 21st june 2017
library(glmnet) # lasso
library(pheatmap) # heatmap
library(e1071) # svm
library(DESeq2) # degs
library(caret) # machine learning methods
library(varSelRF) # variable selection with random forest
library(plyr) # rbind.fill
library(reshape2) # melt
library(tidyverse)
library(leaps) # stepwise feature selection
#-------------------------------------------------------------------------------
#loads bamfiles and ebg data
# 1- bamfiles(Formal class BamFileList) ?
# 2- ebg (GRangesList) is :
#   GRangesList of the human transcriptome

load("data/exons_data.RData")

