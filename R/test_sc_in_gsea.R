

test_sc_in_gsea<-function(targetobj, gs_test=gs_c2_exps, select_gs_names=NULL,padjcutoff=0.01,gsea_cutoff=0.3,negctrlgenelist='NonTargetingControlGuideForHuman'){
  #select_gs_names=grep('^KEGG|^BIOCARTA|^REACTOME|^PID',gs_c2$names,invert = T)
  #padjcutoff=0.01
  if(is.null(select_gs_names)){
    select_gs_names=1:length(gs_test$names)
  }
  r_cs=list()
  r_cs$rpt=list()
  r_cs$dir=c()
  r_cs$pathway=c()

  scalef=getscaledata(targetobj)
  
  for(gs_pw in select_gs_names){
    
    search_ind=which(select_gs_names==gs_pw)
    if(search_ind%%10==1){
      print(paste(search_ind,'/',length(select_gs_names)))
    }
    
    gs_target=gs_test$geneset[gs_pw]
    gs_name=gs_test$names[gs_pw]
    print(gs_name)
    #target_gene='MKI67'
    texp_mat=scalef[rownames(scalef)%in%gs_target[[1]],]
    texp=colSums(texp_mat)
    texp=sort(texp)
    
    #texp_sel_cells=which(names(texp)%in%rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp),'oligo']))
    texp_sel_cells=which(!is.na(targetobj@meta.data[names(texp),'target.gene']))
    texp_withg=texp[texp_sel_cells]
    #texp_guide_ass=bc_dox_uq[names(texp_withg),'oligo']
    #texp_gene_ass=bc_dox_uq[names(texp_withg),'gene']
    texp_guide_ass=targetobj@meta.data[names(texp_withg),'target.gene']
    texp_gene_ass=targetobj@meta.data[names(texp_withg),'geneID']
    
    
    # negative
    
    genes_to_rank=(texp_gene_ass)
    report_neg=get_rank_tables(genes_to_rank,negctrlgenelist=negctrlgenelist)
    
    
    # positive
    
    genes_to_rank=rev(texp_gene_ass)
    report_pos=get_rank_tables(genes_to_rank,negctrlgenelist=negctrlgenelist)
    
    nsel=which(report_neg$gsea>gsea_cutoff & 
                 #report_neg$gene!='TP53'& 
                 (  report_neg$ks_p_adj<padjcutoff | report_neg$ks_negctrl_adj<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'negative'))
      print(report_neg[nsel,])
      r_cs$rpt=rbind(r_cs$rpt,report_neg[nsel,])
      r_cs$dir=c(r_cs$dir,rep('negative',length(nsel)))
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
    nsel=which(report_pos$gsea>gsea_cutoff& 
                 #report_pos$gene!='TP53' & 
                 (report_pos$ks_p_adj<padjcutoff| report_pos$ks_negctrl_adj<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'positive'))
      print(report_pos[nsel,])
      
      
      r_cs$rpt=rbind(r_cs$rpt,report_pos[nsel,])
      r_cs$dir=c(r_cs$dir,rep('positive',length(nsel)))
      
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
  }
  
  report_fall=r_cs$rpt
  report_fall[,'direction']=r_cs$dir
  report_fall[,'pathway']=r_cs$pathway
  report_fall[,'ks_p_adj.2']=p.adjust(report_fall$ks_p,method='fdr')
  report_fall[,'ks_negctrl_adj.2']=p.adjust(report_fall$ks_negctrl,method='fdr')
  return (report_fall)
}
