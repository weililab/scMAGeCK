#convert signature score/pvalue matrix to the format of final results
getsigresult <- function(signature_score, signature_pval) {
  p <- rownames(signature_score)
  output <- data.frame(rep(p, ncol(signature_score)))
  sig_name <- NULL
  for (i in colnames(signature_score)) {
    sig_name <- rbind(sig_name, data.frame(rep(i, nrow(signature_score))))
  }
  sig_sc <- NULL
  for (n in 1:ncol(signature_score)) {
    sig_sc <- rbind(sig_sc, data.frame(signature_score[, n]))
  }
  sig_p <- NULL
  for (n in 1:ncol(signature_pval)) {
    sig_p <- rbind(sig_p, data.frame(signature_pval[, n]))
  }
  output <- cbind(output, sig_name)
  output <- cbind(output, sig_sc)
  output <- cbind(output, sig_p)
  colnames(output) <- paste(c("sgrna", "gene_signature", "LR_score", "p_value"))
  rownames(output) <- (1:nrow(output))
  return(output)
}
TRUE