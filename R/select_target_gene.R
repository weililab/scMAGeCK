select_target_gene<-function(rds_object, bc_frame,perturb_gene, non_target_ctrl, 
                             assay_for_expcor='MAGIC_RNA', # assays used for estimating RNAs
                             min_gene_num=200, max_gene_num=500, logfc.threshold=0.1, min_abs_diff=0.1){
  
  # old method: directly use the "gene" column in Seurat metadata 
  # note: not working for high moi screens
  #deframe=FindMarkers(rds_object,ident.1=perturb_gene,ident.2 = non_target_ctrl,
  #                    group.by = 'gene',logfc.threshold = logfc.threshold)
  if (!'gene'%in% colnames(bc_frame) | ! 'cell'%in% colnames(bc_frame)){
    stop('barcode frame must include columns named gene or cell.')
  }
  cell_perturb=bc_frame[bc_frame$gene%in%perturb_gene,'cell']
  cell_ctrl=bc_frame[bc_frame$gene %in% non_target_ctrl,'cell']
  cell_perturb=cell_perturb[!cell_perturb%in% cell_ctrl]

  cell_perturb=cell_perturb[cell_perturb %in% Cells(rds_object)]
  cell_ctrl=cell_ctrl[cell_ctrl %in% Cells(rds_object)]

  if(length(cell_perturb)<10){
    if(length(cell_perturb)==0){
      stop(paste('No cells express sgRNAs targeting perturbed gene:',perturb_gene,'. Check whether (1) your cell names in barcode file match cell names in your Seurat object; and (2) there are cells that express sgRNAs targeting ',perturb_gene,'.'))
    }else{
      warnings(paste('Only ',length(cell_perturb),'cells express sgRNAs targeting perturbed gene:',perturb_gene))
    }
  }
  if(length(cell_ctrl)<10){
    if(length(cell_ctrl)==0){
      stop(paste('No cells express negative control sgRNAs:',non_target_ctrl,'. Check whether (1) your cell names in barcode file match cell names in your Seurat object; and (2) there are cells that express sgRNAs targeting ',non_target_ctrl,'.'))
    }else{
      warnings(paste('Only ',length(cell_ctrl),'cells express sgRNAs targeting perturbed gene:',non_target_ctrl))
    }
  }
  cells_perturb_meta=rep('Others',length(Cells(rds_object)))
  names(cells_perturb_meta)=Cells(rds_object)
  cells_perturb_meta[cell_perturb]='Perturbed'
  cells_perturb_meta[cells_ctrl]='Ctrl'

  rds_object=AddMetaData(rds_object,cells_perturb_meta,col.name='scmageck_perturbed_marker')
  deframe=FindMarkers(rds_object,ident.1='Perturbed',indent.2='Ctrl',
		      group.by='scmageck_perturbed_marker',
                      logfc.threshold = logfc.threshold)
  
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
                                  assay_for_expcor = assay_for_expcor,
      				  cell1s=cell_perturb,cell_ctrl=cell_ctrl)
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
