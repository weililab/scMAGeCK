pre_processRDS <- function(BARCODE, RDS, normalize = TRUE, scale = TRUE) {
  
  # read cell assignment and libray file ####
  bc_dox = NULL
  if (is.character(BARCODE)) {
    bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
  } else {
    bc_dox = BARCODE
  }
  if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
    stop("cell, barcode, or gene column names not found in barcode file.")
  }
  
  guide_count = table(bc_dox$cell)
  ncnt = table(table(bc_dox$cell))
  message(paste("Total barcode records:", nrow(bc_dox)))
  
  # read Seurat RDS file ####
  if (is.character(RDS)) {
    message(paste("Reading RDS file:", RDS))
    targetobj = readRDS(RDS)
  } else {
    targetobj = RDS
  }
  
  # check if names are consistent
  nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
  if (nmatch == 0) {
    message("Cell names in expression matrix and barcode file do not match. Try to remove possible trailing \"-1\"s...")
    if (length(grep("-\\d$", bc_dox[, 1])) > 0) {
      bc_dox[, 1] = sub("-\\d$", "", bc_dox[, 1])
    }
    nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
    if (nmatch == 0) {
      stop("No cell names match in expression matrix and barcode file.")
    }
  }
  
  #pre-process the barcode file for selection
  # replaced by assign_cell_identity function
  message(paste("Assigning cell identity with highest gRNA expressions (you can change it later by rerunning assign_cell_identity function"))
  targetobj = assign_cell_identity(bc_dox,targetobj,ASSIGNMETHOD='largest')
  #message(paste("Integrating the information of sgRNA into RDS file"))
  #L <- unique(bc_dox$cell)
  #output <- NULL
  #for (i in L) {
  #  a <- nrow(subset(bc_dox, cell == i))
  #  if (a > 1) {
  #    m <- matrix(c(i, NA, NA), ncol = 3, dimnames = list(i, c("cell", "sgrna", "gene")))
  #    output <- rbind(output, m) 
  #  } else {
  #    f <- subset(bc_dox, cell == i)
  #    f <- f[, c("cell", "sgrna", "gene")]
  #    output <- rbind(output, f)
  #  }
  #}
  #row.names(output) <- output$cell #substitute the certain pattern
  #output <- output[, -1]
  
  #merge the grna information into metadata
  #mdata <- targetobj@meta.data
  #mdata <- merge(mdata, output, by = 0, all = TRUE)
  #row.names(mdata) <- mdata$Row.names
  #mdata <- mdata[!is.na(mdata$orig.ident),]
  #mdata <- mdata[, -1]
  #targetobj@meta.data <- mdata
  
  # prepare sgRNA matrix
  message(paste("Creating sgRNA expression matrix"))
  bmatrix <- matrix(0, nrow = length(Cells(targetobj)), ncol = length(unique(bc_dox$gene)))
  rownames(bmatrix) <- Cells(targetobj)
  colnames(bmatrix) <- unique(bc_dox$gene)
  celllist <- rownames(bmatrix)
  genelist <- unique(bc_dox$gene)
  #for (c in cellist) {
  #  for (g in genelist) {
  #    tran <- subset(bc_dox, bc_dox$gene == g & bc_dox$cell == c)
  #    if (nrow(tran) > 0) {
  #      tran <- as.numeric(tran$umi_count)
  #      tran <- sum(tran)
  #      bmatrix[c, g] <- tran
  #    } else {
  #      next
  #    }
  #  }
  #}
  ag_frame=aggregate(umi_count~cell + gene,data=bc_dox,sum)
  ag_frame=ag_frame[ag_frame$cell%in%celllist & ag_frame$gene %in%genelist,]
  bmatrix[cbind(ag_frame$cell,ag_frame$gene)]=ag_frame$umi_count

  bmatrix <- t(bmatrix)
  obj.sgrna <- as.sparse(bmatrix)
  targetobj[['sgRNA']] <- CreateAssayObject(counts = obj.sgrna)
  if ( normalize ) {
    targetobj <- NormalizeData(targetobj, assay = "sgRNA", normalization.method = "CLR")
  }
  if ( scale ) {
    targetobj <- ScaleData(targetobj, assay = "sgRNA")
  }

  
  # prepare sgRNA matrix
  message(paste("Creating sgRNA expression matrix at the guide level"))
  sgmatrix <- matrix(0, nrow = length(Cells(targetobj)), ncol = length(unique(bc_dox$barcode)))
  rownames(sgmatrix) <- Cells(targetobj)
  colnames(sgmatrix) <- unique(bc_dox$barcode)
  celllist <- rownames(sgmatrix)
  barcodelist <- unique(bc_dox$barcode)

  ag_frame=aggregate(umi_count~cell + barcode,data=bc_dox,sum)
  ag_frame=ag_frame[ag_frame$cell %in% celllist & ag_frame$barcode %in% barcodelist,]
  sgmatrix[cbind(ag_frame$cell,ag_frame$barcode)]=ag_frame$umi_count

  sgmatrix <- t(sgmatrix)
  obj.sgrna2 <- as.sparse(sgmatrix)
  targetobj[['sgRNA_guides']] <- CreateAssayObject(counts = obj.sgrna2)
  if ( normalize ) {
    targetobj <- NormalizeData(targetobj, assay = "sgRNA_guides", normalization.method = "CLR")
  }
  if ( scale ) {
    targetobj <- ScaleData(targetobj, assay = "sgRNA_guides")
  }
}
