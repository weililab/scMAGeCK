# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

plot_single_genes<-function(targetob,gene1,targetgene,haslog=TRUE,plotfigure=TRUE,ngctrlgene=c('NonTargetingControlGuideForHuman'),indmatrix=NULL){
  # plot single gene ko effect
  # if indmatrix =NULL, use targetobj@metadata to identify cell groups
  # else, use indmatrix where cells that do not bear gene1 target as controls

  if(is.null(indmatrix)){
    cell_ctrl=rownames(targetob@meta.data)[which(targetob@meta.data$geneID%in% ngctrlgene)]
    cell_gene1=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene1)]
    if(length(cell_gene1)==0){
      gene1=sub('-','_',gene1)
      cell_gene1=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene1)]
    }
  }else{
     cell_ctrl=rownames(indmatrix)[indmatrix[,gene1]==FALSE]
     cell_gene1=rownames(indmatrix)[indmatrix[,gene1]==TRUE]
  }
  #cell_gene2=rownames(targetob@meta.data)[which(targetob@meta.data$geneID==gene2)]
  #mg_geneid=paste(gene1,gene2,sep='_')
  #cell_merged=select_pair_genes[[mg_geneid]]
  scalef=getscaledata(targetobj) 
  t_exp=scalef[targetgene,c(cell_ctrl,cell_gene1)]
  if(min(t_exp)<0){
    t_exp=t_exp-(min(t_exp))+0.1
  }
  t_type=c(rep('NegCtrl',length(cell_ctrl)),rep(gene1,length(cell_gene1)))
  t_type=factor(t_type,levels = c('NegCtrl',gene1))
  ds=data.frame(Expression=t_exp,
                Type=t_type, Cells=c(cell_ctrl,cell_gene1))
  
  if(plotfigure){
    p<-ggplot(ds, aes(x=Type, y=Expression)) + 
      geom_violin()+
      geom_point(position=position_jitter(w=0.1,h=0)) +
      ggtitle(paste(targetgene,'expression'))
    if(haslog){
      p=p+  scale_y_log10()
    }
    #ggtitle(paste(ensemblID,geneID))
    print(p)
  }
  
  return (ds)
}
