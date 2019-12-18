
get_sc_signature <- function(targetobj, gs_test = gs_c2_exps, padjcutoff = 0.01, gsea_cutoff = 0.3, 
    negctrlgenelist = "NonTargetingControlGuideForHuman") {
    # select_gs_names=grep('^KEGG|^BIOCARTA|^REACTOME|^PID',gs_c2$names,invert = T) padjcutoff=0.01
    select_gs_names = gs_test$names
    
    scalef = getscaledata(targetobj)
    
    r_cs = matrix(rep(0, length(select_gs_names) * ncol(scalef)), nrow = length(select_gs_names))
    rownames(r_cs) = select_gs_names
    colnames(r_cs) = colnames(scalef)
    for (gs_pw in 1:length(select_gs_names)) {
        
        # search_ind=which(select_gs_names==gs_pw)
        if (gs_pw%%10 == 1) {
            message(paste(gs_pw, "/", length(select_gs_names)))
        }
        
        gs_target = gs_test$geneset[gs_pw]
        gs_name = gs_test$names[gs_pw]
        # message(gs_name) target_gene='MKI67'
        texp_mat = scalef[rownames(scalef) %in% gs_target[[1]], ]
        texp = colMeans(texp_mat)
        r_cs[gs_pw, ] = texp
        # texp=sort(texp)
    }
    return(r_cs)
    
}

TRUE
