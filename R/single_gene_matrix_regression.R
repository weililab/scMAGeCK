# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

# construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes *
# expressed genes)

single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"), 
    indmatrix = NULL, high_gene_frac = 0.01, selected_genes_list = NULL) {
    # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
    # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
    outlier_threshold = 0.95
    rawf = getscaledata(targetobj, scaled = FALSE)
    select_genes = rownames(rawf)[which(rowSums(as.matrix(rawf) != 0) >= ncol(rawf) * high_gene_frac)]
    if (is.null(selected_genes_list) == FALSE) {
        select_genes = select_genes[select_genes %in% selected_genes_list]
        if (length(select_genes) == 0) {
            stop("No genes left for regression. Check your selected gene list.")
        }
    }
    message(paste("Selected genes:", length(select_genes)))
    # browser()
    
    scalef = getscaledata(targetobj)
    
    if (is.null(indmatrix)) {
        select_cells = rownames(targetobj@meta.data)[which(!is.na(targetobj@meta.data$geneID))]
    } else {
        select_cells = rownames(indmatrix)
        select_cells = select_cells[select_cells %in% colnames(scalef)]
    }
    YmatT = scalef[select_genes, select_cells]
    
    Ymat = as.matrix(t(YmatT))  # (cells * expressed genes)
    if (is.null(indmatrix)) {
        tgf = targetobj@meta.data[select_cells, "geneID"]
        tgf[tgf %in% ngctrlgene] = "NegCtrl"
        tgphenotype = as.factor(tgf)
        Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
        rownames(Xmat) = select_cells
        colnames(Xmat) = levels(tgphenotype)
        Xmat[as.matrix(cbind(1:nrow(Xmat), as.numeric(tgphenotype)))] = 1
        Xmat[, "NegCtrl"] = 1  # set up base line
    } else {
        tgf = colnames(indmatrix)
        tgf[tgf %in% ngctrlgene] = "NegCtrl"
        tgphenotype = as.factor(tgf)
        
        Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
        rownames(Xmat) = select_cells
        colnames(Xmat) = levels(tgphenotype)
        for (cnl in colnames(indmatrix)) {
            cellns = which(indmatrix[, cnl] == TRUE)  #make sure indmatrix 
            if (cnl %in% ngctrlgene) {
                Xmat[cellns, "NegCtrl"] = 1
            } else {
                Xmat[cellns, cnl] = 1
            }
        }
        Xmat[, "NegCtrl"] = 1
        
    }  # end if
    
    # remove outliers
    Ymat_outlier = apply(Ymat, 2, function(X) {
        return(quantile(X, probs = outlier_threshold))
    })
    outlier_mat = t(matrix(rep(Ymat_outlier, nrow(Ymat)), ncol = nrow(Ymat)))
    Ymat_corrected = ifelse(Ymat > outlier_mat, outlier_mat, Ymat)
    Ymat = Ymat_corrected
    
    return(list(Xmat, Ymat))
}
TRUE
