scmageck_lr <- function(BARCODE, RDS, NEGCTRL, SELECT_GENE = NULL, LABEL = NULL, PERMUTATION = NULL,
    SIGNATURE = NULL, SAVEPATH = "./", LAMBDA = 0.01, GENE_FRAC = 0.01, SLOT = 'scale.data') {
  if (!is.null(LABEL)) {
    data_label = LABEL
  } else {
    data_label = "sample1"
  }
  
  if (!is.null(PERMUTATION)) {
    n_permutation = as.integer(PERMUTATION)
  } else {
    n_permutation = 10000
  }
  
  if (!is.null(SIGNATURE)) {
    data_signature=SIGNATURE
    message(paste("run_signature: TRUE"))
  } else {
    data_signature = NULL
    message(paste("run_signature: FALSE"))
  }
  
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
  
  # load neg control guides ####
  ngctrlgenelist = strsplit(NEGCTRL, ",")[[1]]
  message(paste("Neg Ctrl guide:", paste(ngctrlgenelist, collapse = ";")))
  
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
  # bc_dox[,1]=sub('-\\d$','',bc_dox[,1])
  
  # convert to ind_matrix ####
  ind_matrix <- frame2indmatrix(bc_dox, targetobj)
  message(paste("Index matrix dimension:", nrow(ind_matrix), ",", ncol(ind_matrix)))
  
  # try to perform matrix regresson on single genes ####
  mat_for_single_reg = single_gene_matrix_regression(targetobj, selected_genes_list = SELECT_GENE, 
      ngctrlgene = ngctrlgenelist, indmatrix = ind_matrix, high_gene_frac = GENE_FRAC, slot = SLOT)
  Xmat = mat_for_single_reg[[1]]
  
  # Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function
  Ymat = mat_for_single_reg[[2]]
  
  # Optional function
  # Get the results based on gmt file
  if(!is.null(data_signature)){
    gmt <- read_gmt_file(data_signature)
    sig_mat <- getsigmat(Ymat, gmt_file = gmt)
    if (ncol(sig_mat) > 0) {
      Amat_sig_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, sig_mat, lambda = LAMBDA, npermutation = n_permutation)
      sig_score = Amat_sig_lst[[1]]
      sig_pval = Amat_sig_lst[[2]]
      sig_re <- getsigresult(signature_score = sig_score, signature_pval = sig_pval)
      sig_re$Fdr <- p.adjust(sig_re$p_value, method = "fdr")
      write.table(data.frame(sig_re), file = file.path(SAVEPATH, paste(data_label, "_signature.txt", sep = "")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      return(list(sig_result=data.frame(sig_re),
             regression_matrix=list(Xmat=Xmat,Ymat=sig_mat)))
    }
  } else {
    # remove values in Y mat
    Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, Ymat, lambda = LAMBDA, npermutation = n_permutation)
    Amat = Amat_pm_lst[[1]]
    Amat_pval = Amat_pm_lst[[2]]
    if(!is.null(SAVEPATH)){
      write.table(data.frame(Perturbedgene = rownames(Amat), Amat), file = file.path(SAVEPATH, paste(data_label,
         "_score.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(data.frame(Perturbedgene = rownames(Amat), Amat_pval), file = file.path(SAVEPATH, paste(data_label,
         "_score_pval.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(list(score=data.frame(Perturbedgene = rownames(Amat), Amat), 
      pvalue=data.frame(Perturbedgene = rownames(Amat), Amat_pval),
      regression_matrix=list(Xmat=Xmat,Ymat=Ymat)))
  }
}
TRUE
