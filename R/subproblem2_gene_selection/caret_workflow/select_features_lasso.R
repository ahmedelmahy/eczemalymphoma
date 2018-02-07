x_train <-  model.matrix(class~., d_with_class_train) [ , -1]
x_test <- model.matrix(class~., d_with_class_test) [ , -1]

lasso_mod_cv = cv.glmnet(x_train, y_train, alpha =1)
bestlam = lasso_mod_cv$lambda.min
#predict(lasso_mod_cv,newdata = x_test,type = "coefficients", s = bestlam )

p_lasso <- ifelse(predict(lasso_mod_cv ,type ="response", s = bestlam , newx =x_test) <.5, 0,1)

table(y_test,p_lasso)
