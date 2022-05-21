# set working directory
setwd(dir="C:/Users/lucas/Documents/Facc/Analyse Multivariée/Travail Individuel/Multivariate-Analysis-Project/Code")

# clearing the environment and the console
rm(list=ls())
cat("\014")

# import user-defined functions
source(file="source.R")

# 1)
load(file="../Data/Low Dimension Data/Case 1.R")

p <- dim(X1)[2]
m <- 250*p

# random direction & projected data
a <- new_direction(Xj=X1, p=p)
proj_X1 <- X1 %*% a

# medcouple & adjusted outlyingness
MC <- medcouple(proj_X1)

