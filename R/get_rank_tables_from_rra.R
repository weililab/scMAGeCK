

get_rank_tables_from_rra<-function(rankexp,bc_dox_u,rrapath=NULL,pcutoff=0.3,tmpprefix=paste('sample_',runif(1,1,10000),sep=''),negctrlgenelist='NonTargetingControlGuideForHuman',more_rra='',negsel=TRUE,possel=TRUE,keeptmp=FALSE){
  rankexp=rankexp[names(rankexp)%in%rownames(bc_dox_u) & !is.na(bc_dox_u[names(rankexp),'barcode'])]
  if(length(rankexp)<3){
    print('Error: cannot find enough cells.')
    return (NULL)
  }
  rankexp=sort(rankexp)
  
  texp_guide_ass=bc_dox_u[names(rankexp),'barcode']
  texp_gene_ass=bc_dox_u[names(rankexp),'gene']
  texp_guide_ass1=paste(texp_guide_ass,1:length(texp_guide_ass),sep='_r')
  
  rra_oframe=data.frame(guide=texp_guide_ass1,gene=texp_gene_ass,
                        list=rep('list',length(texp_guide_ass)),value=rankexp,
                        prob=rep(1,length(texp_guide_ass)),chosen=rep(1,length(texp_guide_ass)))
  low_file=paste(tmpprefix,'_rra_low.txt',sep='')
  write.table(rra_oframe,file=low_file,row.names = FALSE,quote=FALSE,sep='\t')
  
  
  rra_oframe_h=rra_oframe
  rra_oframe_h[,'value']=-1*rra_oframe_h[,'value']
  rra_oframe_h=rra_oframe_h[order(rra_oframe_h[,'value']),]
  high_file=paste(tmpprefix,'_rra_high.txt',sep='')
  write.table(rra_oframe_h,file=high_file,row.names = FALSE,quote=FALSE,sep='\t')
  ngguidefile=paste(tmpprefix,'_negctrl.txt',sep='')
  
  if(!is.null(negctrlgenelist)){
    ngguidelist=texp_guide_ass1[texp_gene_ass%in%negctrlgenelist]
    write.table(ngguidelist,file=ngguidefile,sep='\t',row.names = FALSE,col.names = FALSE,quote = FALSE)
    ngguidecommand=paste('--control',ngguidefile)
  }else{
    ngguidecommand=''
  }
  if(is.null(rrapath)){
    rracommand='RRA'
  }else{
    rracommand=rrapath
  }
  
  rra_low_out=paste(tmpprefix,'_rra_low.out',sep='')
  rra_c=paste(rracommand,'-i', low_file,
              '-o',rra_low_out,
              ngguidecommand,
              '-p',pcutoff,
              '--max-sgrnapergene-permutation 10000 ',more_rra)
  if(negsel){
  print(rra_c)
  system(rra_c,ignore.stdout = TRUE,ignore.stderr = TRUE)
  }
  
  rra_high_out=paste(tmpprefix,'_rra_high.out',sep='')
  rra_c=paste(rracommand,'-i', high_file,
              '-o',rra_high_out,
              ngguidecommand,
              '-p',pcutoff,
              '--max-sgrnapergene-permutation 10000 ',more_rra)
  if(possel){
  print(rra_c)
  system(rra_c,ignore.stdout = TRUE,ignore.stderr = TRUE)
  }
  
  # merge both
  if(negsel){
    frame_l=read.table(rra_low_out,header = TRUE,as.is = TRUE,row.names = 1,na.strings = '')
  }
  if(possel){
    frame_h=read.table(rra_high_out,header = TRUE,as.is = TRUE,row.names = 1,na.strings = '')
  }
  
  if(negsel & !possel){
    system(paste('rm',low_file,rra_low_out))
    if(!is.null(negctrlgenelist)){
      system(paste('rm',ngguidefile))
    }
    return (frame_l)
  }
  if(!negsel & possel){
    system(paste('rm',high_file,rra_high_out))
    if(!is.null(negctrlgenelist)){
      system(paste('rm',ngguidefile))
    }
    return (frame_h)
  }
  report_f=merge(frame_l,frame_h,by=0,suffixes=c('.low','.high'))
  
  if(!keeptmp){
    system(paste('rm',low_file,high_file,rra_low_out,rra_high_out))
    if(!is.null(negctrlgenelist)){
      system(paste('rm',ngguidefile))
    }
  }
  return (report_f)
}
