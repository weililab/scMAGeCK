# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

# perform matrix decomposition
# get A for Y=XA
# A= (X^TX)^-1X^TY
# Y:  (cells * expressed genes)
# X: design matrix, (cells * KO genes), can also be (cells * double KO genes)
# A: (KO genes * expressed genes)
getsolvedmatrix<-function(Xm,Ym,lambda=0.01){
  #Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
  TMmat_g= (t(Xm) %*% Xm) + lambda * diag(ncol(Xm))
  
  Amat_g= solve(TMmat_g) %*% t(Xm) %*% Ym
  return (Amat_g)
}
