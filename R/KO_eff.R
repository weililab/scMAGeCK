KO_eff <- function(lr_result) {
  data <- as.data.frame(lr_result[1])
  rownames(data) <- sapply(rownames(data), tolower)
  rownames(data) <- capFirst(rownames(data))
  L <- as.data.frame(rownames(data))
  L <- subset(L, L$`rownames(data)` != "Negctrl" & L$`rownames(data)` != "Posctrl")
  L <- L$`rownames(data)`
  gene_score <- NULL
  for (i in L) {
    sc <- data[i, i]
    if (is.null(sc)) {
      next
    } else {
      gene_score <- rbind(gene_score, data.frame("gene" = i, "lrsc" = sc))
    }
  }
  
  data <- as.data.frame(lr_result[2])
  rownames(data) <- sapply(rownames(data), tolower)
  rownames(data) <- capFirst(rownames(data))
  L <- as.data.frame(rownames(data))
  L <- subset(L, L$`rownames(data)` != "Negctrl" & L$`rownames(data)` != "Posctrl")
  L <- L$`rownames(data)`
  p.value <- NULL
  for (i in L) {
    sc <- data[i, i]
    if (is.null(sc)) {
      next
    } else {
      p.value <- rbind(p.value, data.frame("gene" = i, "p.value" = sc))
    }
  }
  
  LR_sc <- merge(gene_score, p.value, by = "gene", all = TRUE)
  LR_sc <- adjust(LR_sc)
  LR_sc$logpadj <- -log(LR_sc$padj)
  g1 <- subset(LR_sc, padj <= 0.05)
  
  p <- ggplot(LR_sc, aes(x=lrsc, y=logpadj)) + geom_point(alpha=4/5) + 
    geom_text(aes(label=ifelse(padj < 0.05, as.character(LR_sc$gene),'')), hjust=0,vjust=0, color = "red", size = 3) + 
    geom_point(data = g1, colour = "red")
  p + xlim(range(c(LR_sc$lrsc, LR_sc$logpadj))) + ylim(range(c(LR_sc$lrsc, LR_sc$logpadj))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "coral1") + 
    geom_hline(yintercept = 3, linetype = "dashed", color = "coral1")
}