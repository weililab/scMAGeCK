scmageck_eff_estimate<-function(rds_object, bc_frame, perturb_gene, non_target_ctrl, 
                                perturb_target_gene = NULL,
                                scale_factor=3, 
                                target_gene_min=200,target_gene_max=500,
                                assay_for_cor='MAGIC_RNA',
                                subset_rds=TRUE,
				scale_score=TRUE,
				perturb_gene_exp_id_list=NULL,
				lambda=0){
  
  if (is.character(rds_object)) {
    message(paste("Reading RDS file:", rds_object))
    rds_object = readRDS(rds_object)
  } 

  if (is.character(bc_frame)) {
    bc_frame = read.table(bc_frame, header = TRUE, as.is = TRUE)
  } 
  # only consider cells that contains perturbed gene or negative controls
  if(DefaultAssay(rds_object)!='RNA'){
    warnings('Your default assay is not RNA. Did you forget to set default assay?')
    DefaultAssay(rds_object)='RNA'
  }
  bc_frame$gene = as.character(bc_frame$gene) # convert to string in cases this is a factor
  bc_subset<-subset(bc_frame,gene%in% c(perturb_gene,non_target_ctrl))
  
  
  rds_used=rds_object
  if(subset_rds){
    rds_used=subset(rds_object,subset = gene %in% c(perturb_gene,non_target_ctrl))
    #rds_used$gene=as.character(rds_used$gene) # convert to character
  }
  target_gene_results_list=list()
  if(is.null(perturb_target_gene)){
    # return the frame of target gene
    message(paste('Target gene list not provided. Search for target gene list..'))
    if(!is.null(perturb_gene_exp_id_list)){
      if(length(perturb_gene_exp_id_list)!=length(perturb_gene)){
	  stop(paste('perturb_gene_exp_id_list parameter should have the same length as perturb_gene.'))
      }
    }
    for(gl_pt_i in 1:length(perturb_gene)){
      gl_pt=perturb_gene[gl_pt_i]
      message(paste('Search for', gl_pt, 'target genes..'))
      if(sum(names(rds_used)==assay_for_cor)==0){ # Assays(rds_used) doesn't work
        stop(paste('Cannot found desired assay', assay_for_cor, 'for estimating correlations.'))
      }

      perturb_gene_exp_id=NULL
      if(!is.null(perturb_gene_exp_id_list)){
	perturb_gene_exp_id=perturb_gene_exp_id_list[gl_pt_i]
      }
      
      perturbed_target_gene_list=select_target_gene(rds_used, bc_subset, gl_pt, non_target_ctrl, 
                                                    min_gene_num = target_gene_min, 
                                                    max_gene_num = target_gene_max,
                                                    assay_for_expcor = assay_for_cor,
                                                    perturb_gene_exp_id=perturb_gene_exp_id) 
      perturb_target_gene=c(perturb_target_gene,perturbed_target_gene_list$target_gene_list)
      target_gene_results_list[[gl_pt]]=perturbed_target_gene_list
    }
    
  }
  # use lr for initial estimation
  message(paste('Using scmageck_lr for initial beta score estimation...'))
  #browser()
  lr_res<-scmageck_lr(bc_subset, rds_used,non_target_ctrl,
                      SELECT_GENE = perturb_target_gene, LABEL='targetgene_lr', PERMUTATION = 100,SLOT='data')
  
  
  lr_score<-lr_res[[1]]
  lr_score_pval<-lr_res[[2]]
  
  lr_xmat=lr_res$regression_matrix$Xmat
  lr_ymat=lr_res$regression_matrix$Ymat
  
  lr_score_fix=lr_score[,-1]
  
  
  
  #
  message('Performing quadratic optimizaiton...')
  
  # full test
  tr_y=as.matrix(lr_ymat)
  tr_x=as.matrix(lr_xmat)
  tr_score=as.matrix(lr_score_fix)
  
  #browser()
  
  tr_x_vec=as.vector(tr_x)
  
  #tr_x_mask=tr_x
  #tr_x_mask[,'NegCtrl']=0 # set neg ctrl column to zero so they don't get udpated
  #tr_x_mask_vector=as.vector(tr_x_mask)
  
  #opres<-optim(tr_x_vec,obj_func,gr=obj_funct_d,method='L-BFGS-B',
  #             lower=rep(0,length(tr_x_vec)),
  #             upper = rep(scale_factor,length(tr_x_vec)),
  #             Y=tr_y,beta_score=tr_score,X_ncol=ncol(tr_x),maskv=tr_x_mask_vector)
  
  #tr_x_update=matrix(opres$par,ncol=ncol(tr_x))
  #rownames(tr_x_update)=rownames(tr_x)
  #colnames(tr_x_update)=colnames(tr_x)
  tr_x_update=scmageck_optim_core(tr_x,tr_y,tr_score,scale_factor=scale_factor, lambda=lambda )
  
  tr_x_update = tr_x_update / scale_factor
  # scale
  if(scale_score == TRUE){
    for(gl_pt in perturb_gene){
      if(gl_pt %in% colnames(tr_x_update)){
        max_score=max(tr_x_update[,gl_pt])
        if(max_score < 0.01){
           warnings(paste('The maximum score for a gene',gl_pt,'is too small (<0.01). No scaling is applied to the corresponding scores.'))
           max_score=1.0
        }
        tr_x_update[,gl_pt] = tr_x_update[,gl_pt] / max_score
      }else{
           warnings(paste(gl_pt,'not in the list of scores to be scaled. Skip..'))
      }
    }
  }
  eff_estimate=tr_x_update[,perturb_gene]
  
  rds_used=AddMetaData(rds_used,eff_estimate,col.name=paste(perturb_gene,'eff',sep='_'))
    
  
  optimization_matrix=list(tr_x=tr_x,tr_y=tr_y,beta_score=tr_score)
  return (list(eff_matrix=tr_x_update[,colnames(tr_x_update)!='NegCtrl'],
               rds=rds_used, 
               optimization_matrix=optimization_matrix, 
               target_gene_search_result=target_gene_results_list))
}
TRUE
