# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

frame2indmatrix<-function(bc_d,targetobj){

  scalef=getscaledata(targetobj)
  colnames(scalef) = sub('-\\d$','',colnames(scalef))
  rnm=unique(bc_d$cell)
  cnm=unique(bc_d$gene)
  rnm=rnm[!is.na(rnm)]
  rnm=rnm[rnm%in%colnames(scalef)]
  cnm=cnm[!is.na(cnm)]
  ind_matrix=matrix(rep(FALSE,length(rnm)*length(cnm)),nrow=length(rnm))
  rownames(ind_matrix)=rnm
  colnames(ind_matrix)=cnm
  for(si in 1:nrow(bc_d)){
    t_r=bc_d[si,'cell']
    t_c=bc_d[si,'gene']
    if((t_r%in%rnm) & (t_c%in% cnm)){
      ind_matrix[t_r,t_c]=TRUE
    }
  }
  return (ind_matrix)
}
