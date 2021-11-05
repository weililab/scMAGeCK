# Convert guide matrix to triplets of (cell, barcode, count)
guidematrix_to_triplet<- function(count_mat, RDS) {


  dgc_mat=as(count_mat,'dgTMatrix')
  cell_names=colnames(dgc_mat)
  guide_names=rownames(dgc_mat)


  if (is.character(RDS)) {
    message(paste("Reading RDS file:", RDS))
    targetobj = readRDS(RDS)
  } else {
    targetobj = RDS
  }
  
  if (sum(cell_names %in% Cells(RDS))==0) {
    stop("Cell names in guide matrix do not match those in Seurat object. Are columns cell names?")
  }
  if (sum(cell_names %in% Cells(RDS)) < 0.01*length(cell_names)) {
    warning("Less than 1% of cells in guide matrix can be found on Seurat object.")
  }

  bc_frame=data.frame(cell=cell_names[dgc_mat@j+1],
                    barcode=guide_names[dgc_mat@i+1],
                    read_count=dgc_mat@x,
                    umi_count=dgc_mat@x)
  return (bc_frame)
}
TRUE
