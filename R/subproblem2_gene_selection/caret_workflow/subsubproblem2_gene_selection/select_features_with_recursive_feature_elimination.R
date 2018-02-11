# using oob error for recurssive feature elimination
rf.vs <- varSelRF(d_train, y_train_factor)
oobrf <- rf.vs$selected.vars
# prediction
gl = glm(formula(paste0("class ~ ",paste0(oobrf,collapse = " + "))), data = d_with_class_train,
              family=binomial)

p_step <- ifelse(predict(gl, newdata = d_test, type="response") < .5,0,1)
table(y_test, p_step)

#         p_step
#           0  1
# y_test 0  23 6
#        1  1  1

