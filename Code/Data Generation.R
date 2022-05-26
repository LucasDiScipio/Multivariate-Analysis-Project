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

# LOW DIMENSION DATA (p < 10)
# ----------------------------

# 1.a)
p <- 2
n_training <- 250 # n_1 = n_2 = 250
n_testing <- n_training/5
Omega <- diag(p)
alpha <- c(5,5)
mu_X1 <- c(0,0)
mu_X2 <- c(2,2)

# training data
X1_training <- rmsn(n_training, xi=mu_X1, Omega, alpha)
X2_training <- rmsn(n_training, xi=mu_X2, Omega, alpha)

# testing data
X1_testing <- rmsn(n_testing, xi=mu_X1, Omega, alpha)
X2_testing <- rmsn(n_testing, xi=mu_X2, Omega, alpha)

# 1.b)
outliers_percentage <- 0.1

# training data
outliers_number <- n_training*outliers_percentage
X1_training_outliers <- X1_training
X1_training_outliers[1:outliers_number,] <- mvrnorm(n=outliers_number, mu=c(-2.5, -2.5), 0.2*diag(p))

# testing data
outliers_number <- n_testing*outliers_percentage
X1_testing_outliers <- X1_testing
X1_testing_outliers[1:outliers_number,] <- mvrnorm(n=outliers_number, mu=c(-2.5, -2.5), 0.2*diag(p))

save(X1_training, X1_testing, X2_training, X2_testing, X1_training_outliers, X1_testing_outliers,
     p, n_training, n_testing, file="../Data/Low Dimension Data/Case 1.R")

