# set working directory
setwd("C:/Users/lucas/Documents/Facc/Analyse Multivariée/Travail Individuel/Multivariate-Analysis-Project/Code")

# clearing the environment and the console
rm(list=ls())
cat("\014")

# Multivariate Normal Distribution
library(MASS)

# Multivariate Skew-Normal Distribution
# install.packages("sn")
library(sn)

# LOW DIMENSION DATA (p <= 10)
# ----------------------------

# 1.a)
p <- 2
n <- 250 # n_1 = n_2 = 250
Omega <- diag(p)
alpha <- c(5,5)
mu_X1 <- c(0,0)
mu_X2 <- c(2,2)
X1 <- rmsn(n, xi=mu_X1, Omega, alpha)
X2 <- rmsn(n, xi=mu_X2, Omega, alpha)

# 1.b)
outliers_percentage = 0.1
outliers_number = n*outliers_percentage
X1_outliers <- X1
X1_outliers[1:outliers_number,] <- mvrnorm(n=outliers_number, mu=c(-2.5, -2.5), 0.2*diag(p))

save(X1, X2, X1_outliers, file="../Data/Low Dimension Data/Case 1.R")