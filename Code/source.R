
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


new_direction <- function(Xj,p){
# ------------------------------------
# GENERATES NEW RANDOM DIRECTION
# 
# INPUTS:
# the dataset (n by p matrix)
# p: the dimension of the dataset
#
# OUTPUTS:
# A p by 1 random direction vector 
# ------------------------------------
    
  # randomly draw p observations from the dataset
  samples <- Xj[sample(nrow(Xj), p),]
  
  # vector(s) between p observations
  vectors <- samples[1,] - samples[2,]
  
  # direction
  X <- Nullspace(matrix(vectors, nrow=1, ncol=p))
  
  # the orthogonal direction is randomly drawn from the basis
  return( matrix(X[,sample(ncol(X), 1)], nrow=p, ncol=1) )
}


h <- function(med_Xj, x_i, x_l) {
# ----------------------------------------
# COMPUTES THE h() KERNEL OF THE MEDCOUPLE
# 
# INPUTS:
# med_Xj: the sample median 
# x_i & x_j: observations from the Xj group 
#
# OUTPUTS:
# the kernel of the medcouple
# ----------------------------------------
  
  return( (x_l-med_Xj - med_Xj-x_i)/(x_l-x_i) )
}


medcouple <- function(Xj){
# ----------------------
# COMPUTES THE MEDCOUPLE
# 
# INPUTS:
# Xj: the dataset (n by p matrix) 
#
# OUTPUTS:
# the medcouple
# ----------------------
  
  med_Xj <- median(Xj)
  
  H <- vector("numeric", 0)
  
  for (i in 1:dim(Xj)[1]){
    
    for (j in 1:dim(Xj)[1]){

      if( Xj[i] < med_Xj & med_Xj < Xj[j] ){
        
        H <- append(H, h(med_Xj, Xj[i], Xj[j]))
        
      } 
      
      else if ( Xj[j] < med_Xj & med_Xj < Xj[i] ) {
        
        H <- append(H, h(med_Xj, Xj[j], Xj[i]))
        
      }
    }
  }
  
  return( median(H) )
}














