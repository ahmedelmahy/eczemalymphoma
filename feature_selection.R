#-------------------------------------------------------------------------------
# feature selection with stepwise function--------------------------------------
# doesn't work that way
fullmod = glm(class ~., data = new.data,family=binomial)

backwards = step(fullmod)




