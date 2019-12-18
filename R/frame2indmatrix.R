# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

frame2indmatrix <- function(bc_d, targetobj) {
  
  rnm = unique(bc_d$cell)
  cnm = unique(bc_d$gene)
  scalef = getscaledata(targetobj)
  message(paste(length(rnm), "..."))
  message(paste(ncol(scalef), "..."))
  # if(sum(rnm%in%colnames(scalef))==0){ message('Cell names in expression matrix and barcode file do
  # not match. Try to remove possible trailing '-1's...')
  # if(length(grep('-\\d$',colnames(scalef)))>0){ colnames(scalef) =
  # sub('-\\d$','',colnames(scalef)) } if(length(grep('-\\d$',rnm))>0){ rnm =
  # sub('-\\d$','',rnm) } }
  rnm = rnm[!is.na(rnm)]
  rnm = rnm[rnm %in% colnames(scalef)]
  if (length(rnm) == 0) {
    stop("Cell names do not match in expression matrix and barcode.")
  }
  cnm = cnm[!is.na(cnm)]
  ind_matrix = matrix(rep(FALSE, length(rnm) * length(cnm)), nrow = length(rnm))
  rownames(ind_matrix) = rnm
  colnames(ind_matrix) = cnm
  row <- bc_d[, 'cell']
  col <- bc_d[, 'gene']
  test <- (row %in% rnm) & (col %in% cnm)
  idx <- cbind(row[test], col[test])
  ind_matrix[idx]  <- TRUE
  return(ind_matrix)
}
TRUE
