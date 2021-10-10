# Get the gene signature expression matrix from Ymat
getsigmat <- function(Ymat, gmt_file) {
  colgmt <- colnames(gmt_file)
  sig_mat <- data.frame(row.names = rownames(Ymat))
  for (num in (1:ncol(gmt_file))) {
    genes <- as.character(gmt_file[, num])
    if (any(genes %in% colnames(Ymat))) {
      genes <- genes[genes %in% colnames(Ymat)]
      Ymat_sig <- as.data.frame(Ymat[, genes])
      Ymat_sig$m <- rowMeans(Ymat_sig)
      sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m)) # identify whether the genome is mouse or human
    } else {
      genes <- capitalize(tolower(genes))
      if(any(genes %in% colnames(Ymat))) {
        genes <- genes[genes %in% colnames(Ymat)]
        Ymat_sig <- as.data.frame(Ymat[, genes])
        Ymat_sig$m <- rowMeans(Ymat_sig)
        sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m))
      } else {
        message(paste(colnames(gmt_file)[num], "can not found in this dataset"))
        colgmt <- subset(colgmt, colgmt != colnames(gmt_file)[num])
        next
      }
    }
  }
  if (ncol(sig_mat) > 0) {
    colnames(sig_mat) <- colgmt
    sig_mat <- as.matrix(sig_mat)
  } else {
    message("No signatures can be found in this dataset")
  }
  return(sig_mat)
}
TRUE
