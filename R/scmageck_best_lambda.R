scmageck_best_lambda <- function(
    rds_object,
    bc_frame,
    non_target_ctrl,
    lambda_seq,
    pseudogene,
    perturb_gene,
    pseudogene_num = 250
) {
  if (is.character(rds_object)) {
    message(paste("Reading RDS file:", rds_object))
    rds_object = readRDS(rds_object)
  } 
  
  rds_object_meta = rds_object@meta.data
  rand_df = rds_object_meta[sample(nrow(rds_object_meta), size=pseudogene_num), ]
  rds_object_meta[rownames(rand_df), "gene"] <- pseudogene
  rds_object_meta[rownames(rand_df), "sgrna"] <- pseudogene
  rds_object <- AddMetaData(rds_object, metadata = rds_object_meta$gene, col.name = "gene")
  rds_object <- AddMetaData(rds_object, metadata = rds_object_meta$sgrna, col.name = "sgrna")
  
  if (is.character(bc_frame)) {
    bc_frame = read.table(bc_frame, header = TRUE, as.is = TRUE)
  } 
  
  bc_frame2 <- bc_frame
  for (x in 1:length(bc_frame$cell)) {
    bc_frame2$gene[x] = ifelse(bc_frame2$cell[x] %in% rownames(rand_df), pseudogene, bc_frame2$gene[x])
  } 
  
  fp_combine = c()
  
  for (x in lambda_seq) {
    print(paste0("lambda = ", x))
    eff_object2 <- scmageck_eff_estimate(rds_object, bc_frame2, perturb_gene,
                                         non_target_ctrl, assay_for_cor = "RNA", lambda = x)
    eff_estimat2=eff_object2$eff_matrix
    rds_subset2=eff_object2$rds
    
    test_meta=rds_subset2@meta.data[,c('gene', paste(perturb_gene, '_eff',sep=''))]
    colnames(test_meta)[2] <- "eff"
    test_meta$scmageck_class.global <- ifelse(test_meta$eff > 0.5, "KO", "NP")
    
    test_meta2=test_meta[test_meta$gene==perturb_gene,]
    output = prop.table(table(test_meta2$scmageck_class.global))
    
    fp = ifelse(names(output[1]) == "KO", as.numeric(output[1]), 0)
    
    fp_combine=c(fp_combine, fp)
  }
  
  names(fp_combine) = lambda_seq
  
  fp_combine_df <- as.data.frame(fp_combine)
  colnames(fp_combine_df)[1] <- "fp"
  fp_combine_df$lambda_seq = log10(lambda_seq)
  colnames(fp_combine_df)[2] <- "log10(lambda)"

  plot(fp_combine_df$`log10(lambda)`, fp_combine_df$fp,
       ylab = "False Positive Rate", xlab = "Log(lambda)")
  abline(h = 0.10, lty=3, col="blue")
  
  return(fp_combine_df)
}
