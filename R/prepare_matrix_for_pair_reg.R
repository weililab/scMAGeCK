# function definitions #####
# version: 02-22-2019
# should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


prepare_matrix_for_pair_reg<-function(targetobj,bc_dox,Xmat,Ymat,Amat,cell_cutoff=4,ngctrlgene=c('NonTargetingControlGuideForHuman')){
  # this function returls (X, Y) for paired KO regression
  # Xmat, Ymat, Amat: these are regression for single genes
  # the gene and cell column in bc_dox is used to determine the target of each cell
  bc_dox_nonuq=bc_dox[!is.na(bc_dox$gene),]
  dupsq=bc_dox_nonuq[duplicated(bc_dox_nonuq$cell),1]
  bc_dox_nonuq=bc_dox_nonuq[bc_dox_nonuq[,1]%in%dupsq,]
  print(paste('non_uq genes:',nrow(bc_dox_nonuq)))
  
  # get the targeting genes per cell
  cell_list=list()
  for(i in 1:nrow(bc_dox_nonuq)){
    cellname=bc_dox_nonuq[i,1]
    targetgene=bc_dox_nonuq[i,'gene']
    if(cellname %in% names(cell_list)){
      cell_list[[cellname]]=append(cell_list[[cellname]],targetgene)
    }else{
      cell_list[[cellname]]=c(targetgene)
    }
  }
  print(paste('finish compiling',length(cell_list),'cells'))
  #browser()
  # filter pairs
  library(utils)
  pair_genes=list()
  
  for(si in 1:length(cell_list)){
    gls=unique(cell_list[[si]])
    gls_cell_name=names(cell_list)[si]
    if(length(gls)<2){
      next
    }
    if(length(gls)>2){
      #print(paste(length(gls),'combinations...'))
    }
    zi_cb=combn(sort(gls),2)
    for(zi in ncol(zi_cb)){
      zi_c=zi_cb[,zi]
      zi_c_id=paste(zi_c,collapse = '_')
      if(zi_c_id %in% names(pair_genes)){
        pair_genes[[zi_c_id]]=append(pair_genes[[zi_c_id]],gls_cell_name)
      }else{
        pair_genes[[zi_c_id]]=c(gls_cell_name)
      }
    }
  }
  print(paste('finished calculating ',length(pair_genes),'pairs'))

 
  scalef=getscaledata(targetobj)
 
  pair_genes_length=unlist(lapply(pair_genes,length))
  
  select_pair_genes=names(pair_genes_length)[pair_genes_length>=cell_cutoff]
  
  # construct Y and X matrix with only double KO
  select_pair_genes=pair_genes[(select_pair_genes)]
  
  cells_combine=unique(unlist(select_pair_genes))
  cells_combine=cells_combine[cells_combine%in% colnames(scalef)]
  
  print(paste('gene pairs left:',length(cells_combine)))
  
  # construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes * expressed genes)
  #select_genes=rownames(targetobj@raw.data)[ which(rowSums(targetobj@raw.data!=0)>ncol(targetobj@raw.data)/100)]
  select_genes=colnames(Ymat)
  select_cells=rownames(Ymat)
  YmatT_db=scalef[select_genes,cells_combine]
  
  Ymat_db=as.matrix(t(YmatT_db)) # (cells * expressed genes)
  #tgphenotype=targetobj@meta.data[select_cells,'geneID']
  #tgf=targetobj@meta.data[select_cells,'geneID']
  #tgf[tgf%in%ngctrlgene]='NegCtrl'
  #tgphenotype=as.factor(tgf)
  tgphenotype=as.factor(colnames(Xmat))
  # need to calculate residue
  Xmat_db_res=matrix(rep(0,length(cells_combine)*length(unique(tgphenotype))),nrow=length(cells_combine))
  
  rownames(Xmat_db_res)=cells_combine
  colnames(Xmat_db_res)=levels(tgphenotype)
  #cells_combine[as.matrix(cbind(1:nrow(cells_combine),as.numeric(tgphenotype)))]=1
  #browser()
  for(i in 1:nrow(bc_dox_nonuq)){
    t_cellname=bc_dox_nonuq[i,1]
    t_gene_target=bc_dox_nonuq[i,'gene']
    if(t_cellname%in%cells_combine & t_gene_target %in% colnames(Xmat_db_res)){
      Xmat_db_res[t_cellname,t_gene_target]=1
    }
  }
  # residule of Xmat * A =Ymat
  Ymat_db_residule=Ymat_db- Xmat_db_res %*% Amat
  
  print('creating matrix for paired gene...')
  # now, create X matrix
  
  Xmat_db=matrix(rep(0,length(cells_combine)*length(names(select_pair_genes))),nrow=length(cells_combine))
  rownames(Xmat_db)=cells_combine
  colnames(Xmat_db)=names(select_pair_genes)
  
  for ( i in 1:length(select_pair_genes)){
    t_pair_name=names(select_pair_genes)[i]
    for(cn in select_pair_genes[[i]]){
      if(cn %in% cells_combine){
        Xmat_db[cn,t_pair_name]=1
      }
    }
  }
  
  return (list(Xmat_db,Ymat_db_residule,select_pair_genes))
}
