sel_lr <- function(GENE, lr_result, CUTOFF = 0.05, ADJ = "fdr"){
  target_gene_list = strsplit(GENE, ",")[[1]]
  target_gene_list = trimws(target_gene_list)
  lr_sc <- as.data.frame(lr_result[1])
  lr_p <- as.data.frame(lr_result[2])
  tem_genelist <- target_gene_list
  for (t in target_gene_list) {
    if(!(t %in% colnames(lr_sc))){
      message(t, " couldn't be found in the dataset")
      tem_genelist <- tem_genelist[tem_genelist != t]
    }
  }
  if (length(tem_genelist >= 1)) {
    sel_sc <- lr_sc[tem_genelist]
    sel_sc <- setNames(sel_sc, paste0(colnames(sel_sc), "_sc"))
    sel_p <- lr_p[tem_genelist]
    sel_p <- setNames(sel_p, paste0(colnames(sel_p), "_p"))
    for (i in colnames(sel_p)) {
      sel_p[, paste0(i,"adj")] <- p.adjust(sel_p[[i]], method = ADJ)
    }
    sel_result <- cbind(sel_sc, sel_p)
    tem_result <- sel_result
    if (length(tem_genelist) == 1) {
      tem_result$FDR <- ifelse(tem_result[[3]] <= CUTOFF, paste("<=", CUTOFF), "other")
      tem_result$sgRNA <- rownames(tem_result)
      tem_result <- tem_result[order(tem_result[[1]], decreasing = FALSE), ]
      tem_result$sgRNA <- factor(tem_result$sgRNA, levels = tem_result$sgRNA)
      ggplot(tem_result, aes(`sgRNA`, y=tem_result[[1]], label=tem_result[[1]])) + 
        geom_bar(stat='identity', aes(fill=FDR), width=.5) + coord_flip() + 
        ggtitle(paste(tem_genelist,"slection plot")) + ylab("LR_score") + theme_bw() 
    } else {
      tem_p <- tem_result[, grep("_padj", colnames(tem_result))]
      tem_p <- tem_p[apply(tem_p, MARGIN = 1, function(x) any(x <= CUTOFF)), ]
      if (nrow(tem_p) == 0) {
        message("Didn't find any sgRNAs that significantly affect the given genelist:")
      } else {
        colnames(tem_p) <- sub("_padj", "", colnames(tem_p))
        fal_result <- NULL
        for (i in colnames(tem_p)) {
          tem <-tem_p[i]
          tem <- subset(tem, tem[1] <= CUTOFF)
          if (nrow(tem) != 0) {
            p <- data.frame(genes = colnames(tem), sgRNA = rownames(tem), padj = tem[[1]])
            fal_result <- rbind(fal_result, p)
          } else {
            next
          }
        }
        tem_s <- tem_result[, grep("_sc", colnames(tem_result))]
        colnames(tem_s) <- sub("_sc", "", colnames(tem_s))
        re <- NULL
        for (i in rownames(tem_s)) {
          tem_1 <- as.matrix(tem_s[i, ])
          tem_1 <- t(tem_1)
          p.1 <- data.frame(genes = rownames(tem_1), sgRNA = colnames(tem_1), lr_sc = tem_1[, 1])
          re <- rbind(re, p.1)
        }
        fal_result <- merge(fal_result, re, by = c(1, 2), all = FALSE)
        fal_result[fal_result == 0] <- "0.001"
        fal_result$padj <- as.numeric(fal_result$padj)
        ggplot(fal_result, aes(x=sgRNA, y=-log10(padj), width = 0.3, fill=genes)) +
          geom_bar(position="dodge", stat = "identity") + theme_bw()
      }
    }
  } else {
    message("Please to check the genes' name")
  }
}