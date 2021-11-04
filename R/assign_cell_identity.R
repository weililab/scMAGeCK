# assign cell identity from barcode files
assign_cell_identity <- function(BARCODE, RDS, ASSIGNMETHOD='unique') {

  bc_dox = NULL
  if (is.character(BARCODE)) {
    bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
  } else {
    bc_dox = BARCODE
  }
  # check barcode file
  if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
    stop("cell, barcode, or gene column names not found in barcode file.")
  }
  
  if (ASSIGNMETHOD == "unique") {
	  
    message("Only assign unique cells")
    bc_dox_dup=bc_dox[duplicated(bc_dox$cell),]
    bc_dox_uniq=bc_dox[!bc_dox$cell%in%bc_dox_dup$cell,]

    cell_ass=rep(NA,nrow(RDS@meta.data))
    names(cell_ass)=rownames(RDS@meta.data)

    cell_ass[as.character(bc_dox_uniq$cell)]=bc_dox_uniq$gene
    cell_ass[as.character(unique(bc_dox_dup$cell))]="doublet"

    RDS<-AddMetaData(object = RDS,metadata = cell_ass,col.name = 'gene')
  } else {
    if (ASSIGNMETHOD == "largest") {
      bc_dox=bc_dox[order(bc_dox[,'umi_count'],decreasing = T),] 
      
      bc_dox_nondup=bc_dox[!duplicated(bc_dox$cell),]
      
      cell_ass=rep(NA,nrow(RDS@meta.data))
      names(cell_ass)=rownames(RDS@meta.data)

      cell_ass[as.character(bc_dox_nondup$cell)]=bc_dox_nondup$gene

      RDS<-AddMetaData(object = RDS,metadata = cell_ass,col.name = 'gene')
    }
  }
  return(RDS)
}
TRUE
