obj_funct_d<-function(X,Y,beta_score,X_ncol=X_ncol,maskv=maskv,lambda=0){
  Xm=matrix(X,ncol=X_ncol)
  diffmat=Xm%*% beta_score - Y 
  #browser()
  sprod=diffmat  %*%  t(beta_score) + lambda
  sprod_v=as.vector(sprod)
  sprod_v*maskv
}
TRUE
