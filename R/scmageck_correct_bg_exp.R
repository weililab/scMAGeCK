scmageck_correct_bg_exp<-function(rds,tr_x,tr_y, 
				  # cluster_column = 'seurat_clusters' ,
				  reductions_used = 'umap',
				  min_cell_cluster=10){
  # correcting background expression in tr_y based on clustering information
  # tr_x is an indication matrix {single-cell (rows) * sgRNAs/genes (columns)}
  # tr_y is a expression matrix {single-cell (rows) * target genes (columns)}
  

  tr_y_dup=tr_y
  # clustering information
  if(is.null(Idents(rds))){
    error('Background exp correction but cannot find clustering assignments from meta data. Either do the clustering first or set up Idents of your cells with clustering results')
  }
  cell_name_clusters=as.character(Idents(rds))
  cell_name=Cells(rds)
  names(cell_name_clusters)=cell_name
  #all_clusters=as.character(rds@meta.data[rownames(tr_x),cluster_column])
  all_clusters=cell_name_clusters[rownames(tr_x)]

  # only use neg ctrl guides
  ng_ctrl_cells=rownames(tr_x)[rowSums(tr_x)==1] 
  ng_ctrl_cells_cluster=all_clusters[ng_ctrl_cells]
  other_cells=rownames(tr_x)[rowSums(tr_x)>1]
  other_cells_cluster=all_clusters[other_cells]
  message(paste('Neg ctrl cells:',length(ng_ctrl_cells)))

  # calculate cluster distances
  rds_red=Reductions(rds,slot= reductions_used)
  if(is.null(rds_red)){
    error(paste('Cannot find reduction',rds_red,'from Seurat object, which will be needed for gene expression correction.'))
  }
  embedding_data=rds_red@cell.embeddings
  
  #browser()

  avg_dist_mat=NULL
  avg_dist_mat_rc=c()
  for(this_cluster in unique(all_clusters)){
    cluster_cell = cell_name[cell_name_clusters==this_cluster]
    ebd_avg=colMeans(embedding_data[cluster_cell,])
    if(is.null(avg_dist_mat)){
      avg_dist_mat=ebd_avg
      avg_dist_mat_rc=this_cluster
    }else{
      avg_dist_mat=rbind(avg_dist_mat,ebd_avg)
      avg_dist_mat_rc=c(avg_dist_mat_rc,this_cluster)
    }
  }
  rownames(avg_dist_mat)=avg_dist_mat_rc

  dist_cluster=as.matrix(dist(avg_dist_mat))

  for(this_cluster in unique(all_clusters)){
    message(paste('correcting cluster',this_cluster,'...'))
    ng_ctrl_cells_c=ng_ctrl_cells[which(ng_ctrl_cells_cluster==this_cluster)]
    other_cells_c=other_cells[which(other_cells_cluster==this_cluster)]

    dist_2_other_clusters=sort(dist_cluster[,this_cluster])
    c_index=2
    while(length(ng_ctrl_cells_c)<min_cell_cluster & c_index <= length(dist_2_other_clusters)){
      nearby_cluster=names(dist_2_other_clusters)[c_index]
      ng_ctrl_cells_c2=ng_ctrl_cells[which(ng_ctrl_cells_cluster==nearby_cluster)]
      c_index=c_index+1
      ng_ctrl_cells_c=c(ng_ctrl_cells_c, ng_ctrl_cells_c2)
      message(paste('Not enough cells in cluster',this_cluster,'; adding cluster',nearby_cluster))
    }

    tr_y_avg = colMeans(tr_y[ng_ctrl_cells_c,])
    tr_y_dup_sub=tr_y_dup[c(ng_ctrl_cells_c,other_cells_c),] 
    
    tr_y_dup_sub= t(t(tr_y_dup_sub) - tr_y_avg)

    tr_y_dup[c(ng_ctrl_cells_c,other_cells_c),]=tr_y_dup_sub

  }
  return (tr_y_dup)


}
