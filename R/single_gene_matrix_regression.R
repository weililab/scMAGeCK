# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

# construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes *
# expressed genes)

single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"), 
    indmatrix = NULL, high_gene_frac = 0.01, selected_genes_list = NULL, slot = 'scale.data') {
    # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
    # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
    outlier_threshold = 0.95
    rawf = getscaledata(targetobj, slot = 'counts')
    scalef = getscaledata(targetobj, slot = slot) # may be scale data or data

    select_genes = rownames(scalef)
    select_genes = select_genes [!is.na(select_genes)]
    if (is.null(selected_genes_list) == FALSE) {
        select_genes = select_genes[select_genes %in% selected_genes_list]
        if (length(select_genes) == 0) {
            stop("No genes left after gene list provided. Check your selected gene list.")
        }
    }

    if(nrow(rawf)>0){
        select_genes = select_genes [ select_genes %in% rownames(rawf) ]
        rawf2 = as.matrix(rawf[select_genes,])
        #browser()
        select_genes = rownames(rawf2)[which(rowSums(rawf2 != 0) >= ncol(rawf2) * high_gene_frac)]
        message(paste('Filter genes whose expression is greater than 0 in raw read count in less than',high_gene_frac,'single-cell populations.' ))
    } else {
        message(paste('Cannot find raw read count in Seurat object. Use scaled data instead, and filter genes whose expression is greater than 0 in less than',high_gene_frac,'single-cell populations.' ))
        scalef2 = as.matrix(scalef [select_genes,])
        select_genes = rownames(scalef2)[which(rowSums(scalef2 >= 0) >= ncol(scalef2) * high_gene_frac)]
    }

    if (length(select_genes) == 0) {
        stop("No genes left for regression. Check your selected gene list.")
    }
    message(paste("Selected genes:", length(select_genes)))
    # browser()
    
    
    if (is.null(indmatrix)) {
        select_cells = rownames(targetobj@meta.data)[which(!is.na(targetobj@meta.data$geneID))]
    } else {
        select_cells = rownames(indmatrix)
        select_cells = select_cells[select_cells %in% colnames(scalef)]
    }
    select_cells = select_cells[!is.na(select_cells) & select_cells %in% colnames(scalef)]
    message(paste("Selected cells:", length(select_cells)))
    YmatT = as.matrix(scalef[select_genes, select_cells])
    
    Ymat = t(YmatT)  # (cells * expressed genes)
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
