scmageck_rra <- function(BARCODE, RDS, GENE, RRAPATH = NULL, LABEL = NULL, NEGCTRL = NULL, KEEPTMP = FALSE, 
    PATHWAY = FALSE, SAVEPATH = "./") {
    if (is.null(RRAPATH)) {
        RRAPATH = system.file("bin", "RRA", package = "scMAGeCK")
    }
    message("Checking RRA...")
    if (!file.exists(RRAPATH)) {
        message("RRA is not does not exist! Please check RRA executable file path")
        return(NULL)
    }
    if (!is.null(LABEL)) {
        data_label = LABEL
    } else {
        data_label = "sample1"
    }
    
    # read cell assignment and libray file ####
    bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
    # check barcode file
    if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
        stop("cell, barcode, or gene column names not found in barcode file.")
    }
    
    keep_tmp = KEEPTMP
    message(paste("keep_tmp:", keep_tmp))
    ispathway = PATHWAY
    message(paste("ispathway:", ispathway))
    
    if (!is.null(NEGCTRL)) {
        negctrl_gene = NEGCTRL
    } else {
        negctrl_gene = NULL
    }
    
    # bc_dox[,1]=sub('-\\d$','',bc_dox[,1])
    
    guide_count = table(bc_dox$cell)
    ncnt = table(table(bc_dox$cell))
    
    # only leave cells with unique guides ####
    
    dupsq = bc_dox[duplicated(bc_dox$cell), 1]
    bc_dox_uq = bc_dox[!bc_dox[, 1] %in% dupsq, ]
    rownames(bc_dox_uq) = bc_dox_uq[, 1]
    
    message(paste("Total barcode records:", nrow(bc_dox)))
    message(paste("Unique barcode records:", nrow(bc_dox_uq)))
    
    # test target genes ####
    target_gene_list = strsplit(GENE, ",")[[1]]
    message(paste("Target gene:", paste(target_gene_list, collapse = ";")))
    
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
            bc_dox_uq[, 1] = sub("-\\d$", "", bc_dox_uq[, 1])
            rownames(bc_dox_uq) = bc_dox_uq[, 1]
        }
        nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
        if (nmatch == 0) {
            stop("No cell names match in expression matrix and barcode file.")
        }
    }
    # run RRA ####
    if ("scale.data" %in% names(attributes(targetobj))) {
        scalef = targetobj@scale.data  # for version 2
    } else {
        scalef = GetAssayData(object = targetobj, slot = "scale.data")
    }
    
    if (ispathway == TRUE) {
        for (target_gene in target_gene_list) {
            if (!target_gene %in% rownames(scalef)) {
                message(paste("Error: gene ", target_gene, " not in expression list."))
                quit()
            }
        }
        texp = colMeans(scalef[target_gene_list, ])
        texp = sort(texp)
        texp_withg = texp[names(texp) %in% rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp), "barcode"])]
        other_table = get_rank_tables_from_rra(texp_withg, bc_dox_uq, tmpprefix = paste("sample_", runif(1, 
            1, 10000), sep = ""), rrapath = RRAPATH, keeptmp = keep_tmp, negctrlgenelist = negctrl_gene)
        if (!is.null(SAVEPATH)) {
            write.table(other_table, file = file.path(SAVEPATH, paste(data_label, "_PATHWAY", "_RRA.txt", 
                sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
        }
        return(other_table)
    } else {
        # treat genes separately
        for (target_gene in target_gene_list) {
            if (!target_gene %in% rownames(scalef)) {
                message(paste("Warning: gene ", target_gene, " not in expression list."))
                next
            } else {
                message(paste("Testing gene ", target_gene, "..."))
            }
            texp = scalef[target_gene, ]
            texp = sort(texp)
            texp_withg = texp[names(texp) %in% rownames(bc_dox_uq) & !is.na(bc_dox_uq[names(texp), "barcode"])]
            other_table = get_rank_tables_from_rra(texp_withg, bc_dox_uq, tmpprefix = paste("sample_", 
                runif(1, 1, 10000), sep = ""), rrapath = RRAPATH, keeptmp = keep_tmp, negctrlgenelist = negctrl_gene)
            if (!is.null(SAVEPATH)) {
                write.table(other_table, file = paste(SAVEPATH, data_label, "_", target_gene, "_RRA.txt", 
                  sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
            }
            return(other_table)
        }
    }
}
TRUE
