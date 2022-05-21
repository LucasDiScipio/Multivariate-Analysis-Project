# set working directory
setwd(dir="C:/Users/lucas/Documents/Facc/Analyse Multivariée/Travail Individuel/Multivariate-Analysis-Project/Code")

# clearing the environment and the console
rm(list=ls())
cat("\014")

# import user-defined functions
source(file="source.R")

# 1)
load(file="../Data/Low Dimension Data/Case 1.R")

n <- dim(X1)[1]
p <- dim(X1)[2]
m <- 250*p

# random direction & projected data
a <- new_direction(X1, p)
proj_X1 <- X1 %*% a

# medcouple & adjusted outlyingness
med_proj_X1 <- median(proj_X1)
MC <- medcouple(proj_X1, med_proj_X1)

AO <- c(rep(0, n));
for (i in 1:dim(proj_X1)[1]) {
  
  AO[i] <- adjusted_outlyingness_univariate(proj_X1[i], proj_X1, med_proj_X1, MC)  

}




