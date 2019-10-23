# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

plot_gi_genes<-function(targetobj,gene1,gene2,targetgene,select_pair_genes,haslog=TRUE,plotfigure=TRUE,ngctrlgene=c('NonTargetingControlGuideForHuman')){
  if(gene1>gene2){
    tx=gene2
    gene2=gene1
    gene1=tx
  }
  cell_ctrl=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID%in% ngctrlgene)]
  cell_gene1=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID==gene1)]
  cell_gene2=rownames(targetobj@meta.data)[which(targetobj@meta.data$geneID==gene2)]
  mg_geneid=paste(gene1,gene2,sep='_')
  cell_merged=select_pair_genes[[mg_geneid]]
  cell_merged=cell_merged[cell_merged%in%rownames(targetobj@meta.data)]
  
  #browser()

  scalef=getscaledata(targetobj)
  t_exp=scalef[targetgene,c(cell_ctrl,cell_gene1,cell_gene2,cell_merged)]
  if(min(t_exp)<0){
    t_exp=t_exp-(min(t_exp))+0.1
  }
  t_type=c(rep('NegCtrl',length(cell_ctrl)),rep(gene1,length(cell_gene1)),rep(gene2,length(cell_gene2)),rep(mg_geneid,length(cell_merged)))
  t_type=factor(t_type,levels = c('NegCtrl',gene1,gene2,mg_geneid))
  ds=data.frame(Expression=t_exp,
                Type=t_type, Cells=c(cell_ctrl,cell_gene1,cell_gene2,cell_merged))
  
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
