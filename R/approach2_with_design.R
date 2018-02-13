cts_ecz <- counts(dds.ecz)
class_ecz <- rep("ecz",dim(cts_ecz)[2])
condition_ecz <- colData(dds.ecz)@listData$condition

cts_lym <- counts(dds.lym)
class_lym <- rep("lym",dim(cts_lym)[2])
condition_lym <- colData(dds.lym)@listData$condition

cts <- cbind(cts_ecz,cts_lym)
cts$class <- as.factor(c(class_ecz,class_lym))
cts$condition <- as.factor(c(condition_ecz,condition_lym))




design <- model.matrix(~0+class+condition,data = cts)
