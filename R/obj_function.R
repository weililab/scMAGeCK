obj_func<-function(X,Y,beta_score,X_ncol=X_ncol,maskv=maskv,lambda=0){
  Xm=matrix(X,ncol=X_ncol)
  #browser()
  sum( (Xm %*% beta_score - Y)^2) * 0.5 + lambda * sum(Xm)
}
TRUE
