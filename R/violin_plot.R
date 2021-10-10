violin_plot <- function(BARCODE, GENE, sgRNA, RDS, CONTROL = NULL) {
  pbmc <- readRDS(RDS)
  barcode <- read.delim(BARCODE)
  target_gene_list = strsplit(GENE, ",")[[1]]
  target_gene_list = trimws(target_gene_list)
  
  # #pre-process the barcode file
  # L <- unique(barcode$cell)
  # output <- NULL
  # for (i in L) {
  #   a <- nrow(subset(barcode, cell == i))
  #   if (a > 1) {
  #     m <- matrix(c(i, NA, NA), ncol = 3, dimnames = list(i, c("cell", "sgrna", "gene")))
  #     output <- rbind(output, m) 
  #   } else {
  #     f <- subset(barcode, cell == i)
  #     f <- f[, c("cell", "sgrna", "gene")]
  #     output <- rbind(output, f)
  #   }
  # }
  # row.names(output) <- sub("-.*", "", output$cell) #substitute the certain pattern
  # output <- output[, -1]
  # 
  # #merge the grna information into metadata
  mdata <- pbmc@meta.data
  # mdata <- merge(mdata, output, by = 0, all = TRUE)
  # row.names(mdata) <- mdata$Row.names
  # mdata <- mdata[!is.na(mdata$orig.ident),]
  # mdata <- mdata[, -1]
  
  if (length(target_gene_list) == 1) {
    data <- FetchData(pbmc, GENE) #the values indicate log(TPM)
    colnames(data)[1] <- paste("genes")
  } else {
    data <- FetchData(pbmc, target_gene_list[1])
    for (i in target_gene_list[2:length(target_gene_list)]) {
      data_1 <- FetchData(pbmc, i)
      data <- cbind(data, data_1)
    }
    data$genes <- rowMeans(data)
    
  }
  if (!is.null(CONTROL)) {
    if (length(CONTROL) > 1) {
      clus <- NULL
      for (c in CONTROL) {
        clus.me <- subset(mdata, seurat_clusters == c)
        clus <- rbind(clus.me, clus)
      }
    } else {
      clus <- subset(mdata, seurat_clusters == CONTROL)
    }
    data <- merge(data, clus, by = 0, all = FALSE)
    rownames(data) <- data$Row.names
    gene.cell <- rownames(subset(clus, gene == sgRNA))
    cell_ko <- subset(data, data[[1]] %in% gene.cell)
  } else {
    gene.cell <- rownames(subset(mdata, gene == sgRNA)) #subset the cells with specific gene knockout
    cell_ko <- subset(data, rownames(data) %in% gene.cell) #find out specific gene expression
    rownames(cell_ko) <- cell_ko$Row.names
  }
  cell_ko$label <- paste(sgRNA, " gene_KO")
  other <- subset(data, !rownames(data) %in% gene.cell)
  other$label <- paste("other")
  gene <- rbind(cell_ko, other)
    
  
  #wilcoxon test > looking up median difference between two group
  wil <- wilcox.test(cell_ko$genes, other$genes, mu=0, alternative = "two.sided", paired = FALSE)
  eq <- paste0("p_value = ", signif(wil$p.value, digits = 3))
  
  #the base of violin plot
  p <- ggplot(gene, aes(x = label, y = genes, fill = label)) + geom_violin() + geom_violin(trim=FALSE) +
    stat_summary(fun.y=mean, geom="point", size=1, color = "black") + labs(y = "Gene expression_lg(TPM)", subtitle = eq)
  p + scale_fill_brewer(palette = "Blues") + theme_classic() + ggtitle(paste(GENE, "gene expression")) + 
    theme(axis.title.y = element_text(margin = margin(r = 10)), axis.title.x=element_blank())
}