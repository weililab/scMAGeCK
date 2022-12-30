get_exp_cor_frame<-function(rds_object,targetgene,ctrlgene,
                            assay_for_expcor='MAGIC_RNA',
			    cell1s=NULL,cell_ctrl=NULL){
  
# get the expression correlation estimation from target gene KO vs negative controls

  # get cell names
  if(is.null(cell1s)){
    cell1s=rownames(rds_object@meta.data)[which(rds_object@meta.data$gene==targetgene)]
  }
  if(is.null(cell_ctrl)){
    cell_ctrl=rownames(rds_object@meta.data)[which(rds_object@meta.data$gene==ctrlgene)]
  }
  
  # get expression data
  exp_all=GetAssayData(rds_object,assay = assay_for_expcor)
  
  exp_1s=t(as.matrix(exp_all[,cell1s]))
  exp_ctrl=t(as.matrix(exp_all[,cell_ctrl]))
  
  # estimate correlation
  if(!targetgene%in%colnames(exp_1s)){
    stop(paste('Cannot find the expression of ',targetgene,'in expression assay.'))
  }
  cor_1s=cor(exp_1s,(exp_1s[,targetgene]))
  
  cor_ctrl=cor(exp_ctrl,(exp_ctrl[,targetgene]))
  
  # set up frame
  corframe=data.frame(gene=rownames(exp_all),cor_target=cor_1s,cor_ctrl=cor_ctrl,
                      #cor_target_nz=cor_1s_nz,cor_ctrl_nz=cor_ctrl_nz,
                      diff=cor_1s-cor_ctrl)
  
  corframe=corframe[order(corframe$diff,decreasing = T),]
  
  corframe_nona=corframe[!is.na(corframe$diff),]
  corframe_nona=corframe_nona[corframe_nona$gene!=targetgene,]
  
  retls=list(cells_target=cell1s,cells_ctrl=cell_ctrl,
             cordata=corframe_nona)
  
  return (retls)
  
}
TRUE
