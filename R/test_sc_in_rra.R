
test_sc_in_rra<-function(targetobj,  gs_test=gs_c2_exps, select_gs_names=NULL,padjcutoff=0.25,gsea_cutoff=0.3,negctrlgenelist='NonTargetingControlGuideForHuman',rra_path=NULL,tmpprefix='sample1'){
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

 
  bc_dx_uq=colnames(scalef)
  bc_guide_ass=targetobj@meta.data[bc_dx_uq,'target.gene']
  bc_gene_ass=targetobj@meta.data[bc_dx_uq,'geneID']
  
  bc_dx_uq=bc_dx_uq[!is.na(bc_gene_ass)]
  bc_guide_ass=bc_guide_ass[!is.na(bc_gene_ass)]
  bc_gene_ass=bc_gene_ass[!is.na(bc_gene_ass)]
  
  bc_dx_ttb=data.frame(cell=bc_dx_uq,barcode=bc_guide_ass,gene=bc_gene_ass)
  rownames(bc_dx_ttb)=bc_dx_uq
  
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
    
    rept_frm=get_rank_tables_from_rra(texp,bc_dx_ttb,rrapath = rra_path,negctrlgenelist = negctrlgenelist,tmpprefix=tmpprefix )
    
    
    nsel=which(#report_neg$gene!='TP53'& 
      (  rept_frm$FDR.low<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'negative'))
      print(rept_frm[nsel,])
      r_cs$rpt=rbind(r_cs$rpt,rept_frm[nsel,])
      r_cs$dir=c(r_cs$dir,rep('negative',length(nsel)))
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
    nsel=which(#report_pos$gene!='TP53' & 
      (rept_frm$FDR.high<padjcutoff))
    if(length(nsel)>0){
      
      print(paste(gs_name,'positive'))
      print(rept_frm[nsel,])
      
      
      r_cs$rpt=rbind(r_cs$rpt,rept_frm[nsel,])
      r_cs$dir=c(r_cs$dir,rep('positive',length(nsel)))
      
      r_cs$pathway=c(r_cs$pathway,rep(gs_name,length(nsel)))
    }
    
  }
  
  report_fall=r_cs$rpt
  report_fall[,'direction']=r_cs$dir
  report_fall[,'pathway']=r_cs$pathway
  #report_fall[,'ks_p_adj.2']=p.adjust(report_fall$ks_p,method='fdr')
  #report_fall[,'ks_negctrl_adj.2']=p.adjust(report_fall$ks_negctrl,method='fdr')
  return (report_fall)
}
