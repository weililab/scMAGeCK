# virtura_facs_functions
get_rank_tables<-function(genes_to_rank,negctrlgenelist='NonTargetingControlGuideForHuman'){
  score_test=seq(-1.0,1.0,length.out = length(genes_to_rank))
  score_tggene=paste(genes_to_rank,1:length(genes_to_rank),sep='_')
  names(score_test)=score_tggene
  
  candidate_gene=unique(genes_to_rank)
  candidate_gene=candidate_gene[!is.na(candidate_gene)]
  gsea_score=rep(1,length(candidate_gene))
  p_vsall=gsea_score
  p_vsnegctrl=gsea_score
  for(tg_i in 1:length(candidate_gene) ){
    
    tggene=candidate_gene[tg_i]
    # print(paste(tggene,'...'))
    ghit=score_tggene[genes_to_rank==tggene]
    
    gs_score=getoriginalgseascore(score_test,ghit)
    z=max(gs_score)
    gsea_score[tg_i]=z
    
    z=ks.test(score_test[!score_tggene%in%ghit],score_test[ghit])
    p_vsall[tg_i]=z$p.value
    
    
    
    z=ks.test(score_test[genes_to_rank%in%negctrlgenelist],score_test[ghit])
    p_vsnegctrl[tg_i]=z$p.value
    
  }
  
  report_f=data.frame(gene=candidate_gene,gsea=gsea_score,ks_p=p_vsall,ks_negctrl=p_vsnegctrl)
  
  report_f[,'ks_p_adj']=p.adjust(report_f$ks_p,method = 'fdr')
  
  report_f[,'ks_negctrl_adj']=p.adjust(report_f$ks_negctrl,method = 'fdr')
  return (report_f)
}