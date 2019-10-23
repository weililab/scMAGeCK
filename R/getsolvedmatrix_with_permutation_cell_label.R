# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


getsolvedmatrix_with_permutation_cell_label<-function(Xm,Ym,lambda=0.01,npermutation=1000){
  Amat_ret=getsolvedmatrix(Xm,Ym,lambda=lambda)
  Amat_ret_higher=matrix(rep(0,ncol(Amat_ret)*nrow(Amat_ret)),nrow = nrow(Amat_ret))
  rownames(Amat_ret_higher)=rownames(Amat_ret)
  colnames(Amat_ret_higher)=colnames(Amat_ret)
  # permute N times
  # randomly shuffle cell labels
  for(npm in 1:npermutation){
    if(npm%%100==0){
      print(paste('Permutation:',npm,'/',npermutation,'...'))
    }
    cells_shu=sample(rownames(Ym),nrow(Ym))
    Xm_s=Xm[cells_shu,]
    Ym_s=Ym # [cells_shu,]
    rownames(Ym_s)=cells_shu
    Amat_random=getsolvedmatrix(Xm_s,Ym_s,lambda=lambda)
    
    Amat_ret_higher=Amat_ret_higher+ (abs(Amat_random)>abs(Amat_ret))*1.0
    #browser()
  }
  Amat_ret_higher=Amat_ret_higher/npermutation
  return (list(Amat_ret,Amat_ret_higher))
}
