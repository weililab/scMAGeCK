scmageck_optim_core<-function(tr_x,tr_y,tr_score,scale_factor=3, lambda=0){
  
  tr_x_vec=as.vector(tr_x)
  
  tr_x_mask=tr_x
  tr_x_mask[,'NegCtrl']=0 # set neg ctrl column to zero so they don't get udpated
  tr_x_mask_vector=as.vector(tr_x_mask)
  
  opres<-optim(tr_x_vec,obj_func,gr=obj_funct_d,method='L-BFGS-B',
               lower=rep(0,length(tr_x_vec)),
               upper = rep(scale_factor,length(tr_x_vec)),
               Y=tr_y,beta_score=tr_score,X_ncol=ncol(tr_x),maskv=tr_x_mask_vector,lambda=lambda)
  
  tr_x_update=matrix(opres$par,ncol=ncol(tr_x))
  rownames(tr_x_update)=rownames(tr_x)
  colnames(tr_x_update)=colnames(tr_x)
  
  return (tr_x_update)
}
TRUE
