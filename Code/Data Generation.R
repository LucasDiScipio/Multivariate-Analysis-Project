# set working directory
setwd("C:/Users/lucas/Documents/Facc/Analyse Multivariée/Travail Individuel/Multivariate-Analysis-Project")

# clearing the environment and the console
rm(list=ls())
cat("\014")

# Multivariate Normal Distribution
library(MASS)

# Multivariate Skew-Normal Distribution
# install.packages("sn")
# library(sn)

# low dimension data (p <= 10)
# ----------------------------

# 1.a)
# X1
p = 2

n = 250 # n_1 = n_2 = 250

Omega_tild = diag(1,p)

alpha = matrix(5, p, 1)

mu = matrix(0, p, 1)

delta = sqrt(1+t(alpha) %*% Omega_tild %*% alpha)
delta = (Omega_tild %*% alpha) / delta[1]

mu_z = delta*sqrt(2/pi)

cov_Z = Omega_tild - mu_z %*% t(mu_z)

phi_p=mvrnorm(n, mu=mu_z, Sigma=cov_Z)

# X2