#-------------------------------------------------------------------------------
# implement stepwise variable selection on the top 50 genes

# full model with all variables
fullmod = glm(class ~ ., data = d_with_class_train, family=binomial)
# forward and backward variable selection
selected_mod = step(fullmod, direction = "both") # selects only three variables
## same result from
# library(MASS)
# step = stepAIC(fullmod, direction = "both", trace = 0)
# This ends up with 5 genes class ~  X100289462 + X1618 + X1673 + X3627 + X6279 +
# X6318
p_step <- ifelse(predict(selected_mod, newdata = d_test, type="response") < .5,0,1)
table(y_test, p_step)
# with top 17 lym and top 17 ecz
#           p_step
#            0   1
# y_test 0   24  5
#        1   2   0

