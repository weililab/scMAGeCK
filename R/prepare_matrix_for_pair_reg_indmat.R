# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

prepare_matrix_for_pair_reg_indmat <- function(targetobj, ind_matrix, Xmat, Ymat, Amat, cell_cutoff = 4, 
    ngctrlgene = c("NonTargetingControlGuideForHuman")) {
    # this function returns (X, Y) for paired KO regression Xmat, Ymat, Amat: these are regression for
    # single genes the gene and cell column in ind_matrix is used to determine the target of each cell
    
    scalef = getscaledata(targetobj)
    # new code
    s_pair_x = c()
    s_pair_y = c()
    cells_s = rep(0, nrow(ind_matrix))
    names(cells_s) = rownames(ind_matrix)
    for (i in 1:(ncol(ind_matrix) - 1)) {
        if (i%%10 == 1) {
            message(paste(i, "..."))
        }
        for (j in (i + 1):ncol(ind_matrix)) {
            nct = sum(ind_matrix[, i] & ind_matrix[, j])
            if (nct >= cell_cutoff) {
                s_pair_x = append(s_pair_x, i)
                s_pair_y = append(s_pair_y, j)
                cells_s = cells_s + (ind_matrix[, i] & ind_matrix[, j])
            }
        }
    }
    cells_combine = names(cells_s)[which(cells_s > 0)]
    cells_combine = cells_combine[cells_combine %in% colnames(scalef)]
    
    
    # construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes *
    # expressed genes) select_genes=rownames(targetobj@raw.data)[
    # which(rowSums(targetobj@raw.data!=0)>ncol(targetobj@raw.data)/100)]
    select_genes = colnames(Ymat)
    select_cells = rownames(Ymat)
    YmatT_db = scalef[select_genes, cells_combine]
    
    Ymat_db = as.matrix(t(YmatT_db))  # (cells * expressed genes)
    # tgphenotype=targetobj@meta.data[select_cells,'geneID']
    # tgf=targetobj@meta.data[select_cells,'geneID'] tgf[tgf%in%ngctrlgene]='NegCtrl'
    # tgphenotype=as.factor(tgf)
    tgphenotype = as.factor(colnames(Xmat))
    # need to calculate residue
    # Xmat_db_res=matrix(rep(0,length(cells_combine)*length(unique(tgphenotype))),nrow=length(cells_combine))
    # rownames(Xmat_db_res)=cells_combine colnames(Xmat_db_res)=levels(tgphenotype)
    Xmat_db_res = Xmat[cells_combine, ]
    
    # cells_combine[as.matrix(cbind(1:nrow(cells_combine),as.numeric(tgphenotype)))]=1 browser()
    
    # residule of Xmat * A =Ymat
    Ymat_db_residule = Ymat_db - Xmat_db_res %*% Amat
    
    message("creating matrix for paired gene...")
    # now, create X matrix
    
    Xmat_db = matrix(rep(0, length(cells_combine) * length(s_pair_x)), nrow = length(cells_combine))
    rownames(Xmat_db) = cells_combine
    colnames(Xmat_db) = paste(colnames(ind_matrix)[s_pair_x], "_", colnames(ind_matrix)[s_pair_y], sep = "")
    
    message(paste("selected gene pairs:", length(s_pair_x)))
    select_pair_genes = list()
    for (i in 1:length(s_pair_x)) {
        nct = rownames(ind_matrix)[which(ind_matrix[, s_pair_x[i]] & ind_matrix[, s_pair_y[i]])]
        cn = nct[nct %in% cells_combine]
        Xmat_db[cn, i] = 1
        namex = colnames(ind_matrix)[s_pair_x[i]]
        namey = colnames(ind_matrix)[s_pair_y[i]]
        genepairname = paste(namex, "_", namey, sep = "")
        select_pair_genes[[genepairname]] = cn
    }
    
    return(list(Xmat_db, Ymat_db_residule, select_pair_genes))
}
TRUE
