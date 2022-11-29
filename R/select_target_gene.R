select_target_gene<-function(rds_object, perturb_gene, non_target_ctrl, 
                             assay_for_expcor='MAGIC_RNA', # assays used for estimating RNAs
                             min_gene_num=200, max_gene_num=500, logfc.threshold=0.1, min_abs_diff=0.1){
  
  deframe=FindMarkers(rds_object,ident.1=perturb_gene,ident.2 = non_target_ctrl,
                      group.by = 'gene',logfc.threshold = logfc.threshold)
  
  target_gene_list=rownames(deframe)
  corfr_sel=NULL
  message(paste('DE analysis identified ',nrow(deframe),'DE genes that will be used for target genes.'))
  if(nrow(deframe)>max_gene_num){
    deframe=deframe[order(deframe$p_val),]
    target_gene_list=rownames(deframe)[1:max_gene_num]
    message(paste('Select',max_gene_num,'genes'))
  }else{
    if(nrow(deframe)<min_gene_num){
      # select genes from expression correlation
      message('Not enough genes from differential expression. Select genes from correlation...')
      targetls<-get_exp_cor_frame(rds_object,perturb_gene,non_target_ctrl,
                                  assay_for_expcor = assay_for_expcor)
      corfr<-targetls$cordata
      
      corfr_sel=subset(corfr, abs(cor_target)-abs(cor_ctrl)>min_abs_diff)
      corfr_sel=corfr_sel[order(abs(corfr_sel$cor_ctrl),decreasing = T),]
      
      remain_gene=min(nrow(corfr_sel),min_gene_num-nrow(deframe))
      message(paste('Remaining',remain_gene,'genes selected..'))
      target_gene_list=c(target_gene_list,corfr_sel$gene[1:remain_gene])
    }
  }
  
  return (list(target_gene_list=target_gene_list,
               deframe=deframe,
               corframe=corfr_sel))
  
  
}
TRUE
