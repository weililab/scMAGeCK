selectPlot <- function(GENE = NULL, lr_result = NULL, CUTOFF = 0.05, ADJ = "fdr", 
                       RRA_re1 = NULL, RRA_re2 = NULL, TYPE = select.type, QUALITY = 10) {
  if (TYPE == "lr") {
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
      temresult <- sel_result
      if (length(tem_genelist) == 1) {
        temresult$FDR <- ifelse(temresult[[3]] <= CUTOFF, paste("<=", CUTOFF), "other")
        temresult$sgRNA <- rownames(temresult)
        temresult <- temresult[order(temresult[[1]], decreasing = FALSE), ]
        temresult$sgRNA <- factor(temresult$sgRNA, levels = temresult$sgRNA)
        ggplot(temresult, aes(`sgRNA`, y=temresult[[1]], label=temresult[[1]])) + 
          geom_bar(stat='identity', aes(fill=FDR), width=.5) + coord_flip() + 
          ggtitle(paste(tem_genelist,"slection plot")) + ylab("LR_score") + theme_bw() 
      } else {
        tem_p <- temresult[, grep("_padj", colnames(temresult))]
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
          tem_s <- temresult[, grep("_sc", colnames(temresult))]
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
  } else {
    if (TYPE == "rra") {
    na <- paste("<=", CUTOFF, sep = "")
    pos_otx2 <- subset(RRA_re1, p.high < p.low)
    pos_otx2$FDR.n <- ifelse(pos_otx2$FDR.high <= CUTOFF & pos_otx2$goodsgrna.high >= CUTOFF, paste(na), "other")
    pos_otx2$sel_score1 <- -log(pos_otx2$p.high)
    pos_otx2 <- pos_otx2[c("Row.names", "sel_score1", "FDR.n")]
    colnames(pos_otx2)[1] <- paste("markers")
    
    neg_otx2 <- subset(RRA_re1, p.low < p.high)
    neg_otx2$FDR.n <- ifelse(neg_otx2$FDR.low <= CUTOFF & neg_otx2$goodsgrna.low >= CUTOFF, paste(na), "other")
    neg_otx2$sel_score1 <- log(neg_otx2$p.low)
    neg_otx2 <- neg_otx2[c("Row.names", "sel_score1", "FDR.n")]
    colnames(neg_otx2)[1] <- paste("markers")
    n_otx2 <- rbind(pos_otx2, neg_otx2)
    
    if (!is.null(RRA_re2)) {
      pos_otx2 <- subset(RRA_re2, p.high < p.low)
      pos_otx2$FDR.p <- ifelse(pos_otx2$FDR.high <= CUTOFF & pos_otx2$goodsgrna.high >= QUALITY, paste(na), "other")
      pos_otx2$sel_score2 <- -log(pos_otx2$p.high)
      pos_otx2 <- pos_otx2[c("Row.names", "sel_score2", "FDR.p")]
      colnames(pos_otx2)[1] <- paste("markers")
      
      neg_otx2 <- subset(RRA_re2, p.low < p.high)
      neg_otx2$FDR.p <- ifelse(neg_otx2$FDR.low <= CUTOFF & neg_otx2$goodsgrna.low >= QUALITY, paste(na), "other")
      neg_otx2$sel_score2 <- log(neg_otx2$p.low)
      neg_otx2 <- neg_otx2[c("Row.names", "sel_score2", "FDR.p")]
      colnames(neg_otx2)[1] <- paste("markers")
      p_otx2 <- rbind(pos_otx2, neg_otx2)
      
      neg_Nanog <- merge(n_otx2, p_otx2, by = "markers")
      rownames(neg_Nanog) <- neg_Nanog$markers
      neg_Nanog$index <- ifelse(neg_Nanog$FDR.p == na | neg_Nanog$FDR.n == na, rownames(neg_Nanog), "other")
      
      ggplot(neg_Nanog, aes(x=sel_score1, y=sel_score2)) + geom_point(aes(color=index)) + labs(factor = element_blank()) + 
        xlim(range(c(neg_Nanog$sel_score1, neg_Nanog$sel_score2))) + ylim(range(c(neg_Nanog$sel_score1, neg_Nanog$sel_score2))) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "coral1") + geom_hline(yintercept = 0, linetype = "dashed", color = "coral1") + 
        geom_text(aes(label=ifelse(neg_Nanog$FDR.p == na | neg_Nanog$FDR.n == na, as.character(rownames(neg_Nanog)),'')), hjust=0.6,vjust=0, size = 3)
      
    } else {
      rownames(n_otx2) <- n_otx2$markers
      colnames(n_otx2)[3] <- "FDR"
      n_otx2 <- n_otx2[order(n_otx2[[2]], decreasing = FALSE), ]
      n_otx2$markers <- factor(n_otx2$markers, levels = n_otx2$markers)
      ggplot(n_otx2, aes(`markers`, y=sel_score1, label=sel_score1)) + geom_bar(stat='identity', aes(fill=FDR), width=.5) + 
        coord_flip()
    }
  }
  }
}