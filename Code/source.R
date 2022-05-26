
Nullspace <- function (A) {
# --------------------------------------------
# FINDS THE NULLSPACE ASSOCIATED WITH MATRIX A
# 
# INPUTS:
# A: m by n matrix
#
# OUTPUTS:
# X: the basis vector(s) of the nullspace
# --------------------------------------------

  m <- dim(A)[1]; n <- dim(A)[2]
  ## QR factorization and rank detection
  QR <- base::qr.default(A)
  r <- QR$rank
  if ((r < min(m, n)) || (m < n)) {
    R <- QR$qr[1:r, , drop = FALSE]
    P <- QR$pivot
    F <- R[, (r + 1):n, drop = FALSE]
    I <- base::diag(1, n - r)
    B <- -1.0 * base::backsolve(R, F, r)
    Y <- base::rbind(B, I)
    X <- Y[base::order(P), , drop = FALSE]
    return(X)
  }
  ## case 1
  return(base::matrix(0, n, 1))
}


new_direction <- function(X, p){
# --------------------------------
# GENERATES NEW RANDOM DIRECTION
# 
# INPUTS:
# X: the dataset (n by p matrix)
# p: the dimension of the dataset
#
# OUTPUTS:
# A p by 1 random direction vector 
# --------------------------------
    
  # randomly draw p observations from the dataset
  samples <- X[sample(nrow(X), p),]
  
  # vector(s) between p observations
  vectors <- samples[1,] - samples[2,]
  
  # direction
  X <- Nullspace(matrix(vectors, nrow=1, ncol=p))
  
  # the orthogonal direction is randomly drawn from the basis
  return( matrix(X[,sample(ncol(X), 1)], nrow=p, ncol=1) )
}


h <- function(med, x_i, x_l) {
# ----------------------------------------
# COMPUTES THE h() KERNEL OF THE MEDCOUPLE
# 
# INPUTS:
# med: the median of the dataset 
# x_i & x_l: observations from the dataset 
#
# OUTPUTS:
# the kernel of the medcouple
# ----------------------------------------
  
  return( (x_l-med - med-x_i)/(x_l-x_i) )
}


medcouple <- function(X, med, n){
# ----------------------------------------
# COMPUTES THE MEDCOUPLE
# 
# INPUTS:
# X: the dataset
# med: the median of the dataset
# n: the number of elements of the dataset
#
# OUTPUTS:
# the medcouple
# ----------------------------------------
  
  H <- vector("numeric", 0)
  
  for (i in 1:n){
    
    for (j in 1:n){

      if( X[i] < med & med < X[j] ){
        
        H <- append(H, h(med, X[i], X[j]))
        
      } 
      
      else if ( X[j] < med & med < X[i] ) {
        
        H <- append(H, h(med, X[j], X[i]))
        
      }
    }
  }
  
  return( median(H) )
}


adjusted_outlyingness_univariate <- function(x, X, med, MC){
# -----------------------------------------------------------
# COMPUTES THE UNIVARIATE ADJUSTED OUTLYINGNESS 
# 
# INPUTS:
# x: an observation from the dataset
# X: the dataset
# med: the median of the dataset
# MC: the medcouple of the dataset
#
# OUTPUTS:
# the univariate adjusted outlyingness of x with respect to X
# -----------------------------------------------------------
  
  # If MC < 0, the AO is computed on the inverted dataset
  if (MC < 0) X <- -X
  
  IQR <- IQR(X)
  
  if (x >= med) {

    Q_3 <- unname(quantile(X, probs=.75))
    c_2 <- max(subset(X, X < Q_3+1.5*exp(3*MC)*IQR))
    return( (x-med)/(c_2-med) )
    
  } else {

    Q_1 <- unname(quantile(X, probs=.25))          
    c_1 <- min(subset(X, X > Q_1-1.5*exp(-4*MC)*IQR))
    return( (med-x)/(med-c_1) )
    
  }
  
}


adjusted_outlyingness_multivariate <- function(X, m){
# -----------------------------------------------------------------
# COMPUTES THE MULTIVARIATE ADJUSTED OUTLYINGNESS 
# 
# INPUTS:
# X: the dataset
# m: the number of directions to be considered
#
# OUTPUTS:
# the multivariate adjusted outlyingness of the observations from X
# -----------------------------------------------------------------

  n <- dim(X)[1]
  AO <- matrix(0, m, n)
  
  for (i in 1:m){
    
    # random directions & projected data
    a <- new_direction(X, p)
    proj_X <- X %*% a
    
    # medcouple & adjusted outlyingness
    med_proj_X <- median(proj_X)
    MC <- medcouple(proj_X, med_proj_X, n)
    for (j in 1:n) AO[i,j] <- adjusted_outlyingness_univariate(proj_X[j], proj_X, med_proj_X, MC)

    print(i)
  }
  
  # sup of the m AOs for each observation
  max_AO <- matrix(0, nrow=n, ncol=1)
  for (i in 1:n)  max_AO[i] <- max(AO[,i])
    
  return(max_AO)
}


outlier_score <- function(AO, n){
# ----------------------------------
# COMPUTES THE OUTLIER SCORES
# 
# INPUTS:
# AO: the set of AOs
# n: the number of AOs
#
# OUTPUTS:
# OS: the outlier scores for each AO
# ----------------------------------

  med <- median(AO)
  MC <- medcouple(AO, med, n)
  
  OS <- matrix(0, nrow=n, ncol=1)
  for(i in 1:n)  OS[i] <- adjusted_outlyingness_univariate(AO[i], AO, med, MC)
  
  return(OS)
}
















