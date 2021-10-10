sel_rra <- function(RRA_re1, RRA_re2 = NULL) {
  n_otx2 <- read.delim(RRA_re1)
  pos_otx2 <- subset(n_otx2, p.high < p.low)
  pos_otx2$FDR.n <- ifelse(pos_otx2$FDR.high <= 0.05 & pos_otx2$goodsgrna.high >=10, "good", "bad")
  pos_otx2$nscore <- -log(pos_otx2$p.high)
  pos_otx2 <- pos_otx2[c("Row.names", "nscore", "FDR.n")]
  colnames(pos_otx2)[1] <- paste("markers")
  
  neg_otx2 <- subset(n_otx2, p.low < p.high)
  neg_otx2$FDR.n <- ifelse(neg_otx2$FDR.low <= 0.05 & neg_otx2$goodsgrna.low >=10, "good", "bad")
  neg_otx2$nscore <- log(neg_otx2$p.low)
  neg_otx2 <- neg_otx2[c("Row.names", "nscore", "FDR.n")]
  colnames(neg_otx2)[1] <- paste("markers")
  n_otx2 <- rbind(pos_otx2, neg_otx2)
  
  if (!is.null(RRA_re2)) {
    p_otx2 <- read.delim(RRA_re2)
    pos_otx2 <- subset(p_otx2, p.high < p.low)
    pos_otx2$FDR.p <- ifelse(pos_otx2$FDR.high <= 0.05 & pos_otx2$goodsgrna.high >=10, "good", "bad")
    pos_otx2$pscore <- -log(pos_otx2$p.high)
    pos_otx2 <- pos_otx2[c("Row.names", "pscore", "FDR.p")]
    colnames(pos_otx2)[1] <- paste("markers")
    
    neg_otx2 <- subset(p_otx2, p.low < p.high)
    neg_otx2$FDR.p <- ifelse(neg_otx2$FDR.low <= 0.05 & neg_otx2$goodsgrna.low >=10, "good", "bad")
    neg_otx2$pscore <- log(neg_otx2$p.low)
    neg_otx2 <- neg_otx2[c("Row.names", "pscore", "FDR.p")]
    colnames(neg_otx2)[1] <- paste("markers")
    p_otx2 <- rbind(pos_otx2, neg_otx2)
    
    neg_Nanog <- merge(n_otx2, p_otx2, by = "markers")
    rownames(neg_Nanog) <- neg_Nanog$markers
    neg_Nanog$index <- ifelse(neg_Nanog$FDR.p == "good" | neg_Nanog$FDR.n == "good", rownames(neg_Nanog), "other")
    
    ggplot(neg_Nanog, aes(x=nscore, y=pscore)) + geom_point(aes(color=index)) + labs(factor = element_blank()) + 
      xlim(range(c(neg_Nanog$nscore, neg_Nanog$pscore))) + ylim(range(c(neg_Nanog$nscore, neg_Nanog$pscore))) + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "coral1") + geom_hline(yintercept = 0, linetype = "dashed", color = "coral1") + 
      geom_text(aes(label=ifelse(neg_Nanog$FDR.p == "good" | neg_Nanog$FDR.n == "good", as.character(rownames(neg_Nanog)),'')), hjust=0.6,vjust=0, size = 3)
    
  } else {
    rownames(n_otx2) <- n_otx2$markers
    colnames(n_otx2)[3] <- "FDR"
    n_otx2 <- n_otx2[order(n_otx2[[2]], decreasing = FALSE), ]
    n_otx2$markers <- factor(n_otx2$markers, levels = n_otx2$markers)
    ggplot(n_otx2, aes(`markers`, y=nscore, label=nscore)) + geom_bar(stat='identity', aes(fill=FDR), width=.5) + 
      coord_flip()
  }
}