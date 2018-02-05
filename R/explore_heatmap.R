source("R/shiny_pheatmap.R")
df_test = new.data[new.data$class == "lym",]
df_train = new.data[new.data$class == "ecz",]


shiny_pheatmap(datalist = list("train" = df_train,"test" = df_test))
