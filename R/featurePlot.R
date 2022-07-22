featurePlot <- function(RDS, TYPE = plot.type, BARCODE = NULL, sgRNA = NULL, GENE = NULL, CONTROL = NULL, GROUP2=NULL, SLOT='data', 
    palette = NULL, label.size = 3, axis.size = 12, title.size = 15, legend.text = 10, fill = "#56B4E9")  {
  
  p<-NULL
  if (TYPE == "Dis") {
    if (is.null(BARCODE)) {
      stop("Please provide the barcode file from previous step")
    }
    # read cell assignment and libray file ####
    bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
    if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
      stop("cell, barcode, or gene column names not found in barcode file.")
    }
    
    # check if names are consistent
    nmatch = sum(bc_dox[, 1] %in% colnames(x = RDS))
    if (nmatch == 0) {
      message("Cell names in expression matrix and barcode file do not match. Try to remove possible trailing \"-1\"s...")
      if (length(grep("-\\d$", bc_dox[, 1])) > 0) {
        bc_dox[, 1] = sub("-\\d$", "", bc_dox[, 1])
      }
      nmatch = sum(bc_dox[, 1] %in% colnames(x = RDS))
      if (nmatch == 0) {
        stop("No cell names match in expression matrix and barcode file.")
      }
    }
    
    L <- unique(bc_dox$cell)
    output <- NULL
    for (i in L) {
      a <- nrow(subset(bc_dox, cell == i))
      m <- matrix(c(i, a), ncol = 2, dimnames = list(i, c("cell", "number_gRNA")))
      output <- rbind(output, m)
    }
  
    output <- as.data.frame(output)
    output$number_gRNA <- as.numeric(output$number_gRNA)
    med <- as.data.frame(colnames(RDS))
    colnames(med) <- paste("cell")
    med <- merge(med, output, by = 1, all = FALSE)
  
    N <- unique(med$number_gRNA)
    output_1 <- NULL
    for (t in N) {
      b <- nrow(subset(med, number_gRNA == t))
      c <- data.frame(number_gRNA = t, ncell = b)
      output_1 <- rbind(output_1, c)
    }
    
    cell0 <- data.frame(number_gRNA = 0, ncell = (ncol(RDS) - nrow(med)))
    output_1 <- rbind(cell0, output_1)
    
    num <- nrow(output) - nrow(med)
    message(paste(num, "cells entered sgRNAs have been filtered out in RDS file"))
    p <- ggplot(output_1, aes(number_gRNA, ncell)) +
      geom_col(width = 0.4, fill = fill, colour = "black") + 
      ggtitle("sgRNA distribution") + geom_text(aes(label = ncell), size = label.size, hjust = 0.5, vjust = .01) + 
      theme_bw() + theme(axis.text = element_text(size = axis.size)) + 
      theme(plot.title = element_text(size = title.size)) + theme(axis.title = element_text(size = title.size))
    
  } else {
    
    if (TYPE == "Vln" | TYPE == "ECDF") {
      
      if (is.null(GENE)) {
        stop("Target gene is missing")
      }
      if (is.null(sgRNA)) {
        stop("sgRNA is missing")
      }
      
     if (is.vector(GENE)) {
       target_gene_list=GENE
     }else{
       target_gene_list = strsplit(GENE, ",")[[1]]
       target_gene_list = trimws(target_gene_list)
     }
      mdata <- RDS@meta.data

      if (sum(colnames(mdata)=='gene') == 0) {
	  stop("Please assign single cells with gene identity in metadata first.")
      }
      
      if (length(target_gene_list) == 1) {
        data <- FetchData(RDS, GENE, slot = SLOT)  #the values indicate log(TPM)
        colnames(data)[1] <- paste("genes")
      } else {
        data <- FetchData(RDS, target_gene_list[1], slot = SLOT)
        for (i in target_gene_list[2:length(target_gene_list)]) {
          if (i %in% rownames(RDS) ) {
            data_1 <- FetchData(RDS, i, slot = SLOT)
            data <- cbind(data, data_1)
          }else{
            message(paste('Warning:',i,'is missing in expression assays. Skip this gene.'))
          }
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
      if( is.null(GROUP2) ) {
        other <- subset(data, !rownames(data) %in% gene.cell)
        other$label <- paste("other")
      } else {
        other <- subset(data, rownames(data) %in% rownames(subset(mdata, gene == GROUP2)))
        other$label <- paste(GROUP2)
      }
      gene <- rbind(cell_ko, other)
      
      
      #wilcoxon test > looking up median difference between two group
      wil <- wilcox.test(cell_ko$genes, other$genes, mu=0, alternative = "two.sided", paired = FALSE)
      eq <- paste0("p_value = ", signif(wil$p.value, digits = 3))
      #the base of violin plot
      if (TYPE == "Vln"){
          p <- ggplot(gene, aes(x = label, y = genes, fill = label)) + geom_violin() + 
            # geom_violin(trim=FALSE) +  ## trim the figure
            stat_summary(fun.y=mean, geom="point", size=1, color = "black") + labs(y = "Gene expression_lg(TPM)", subtitle = eq)
            
          if (!is.null(palette)) {
            p <- p + scale_fill_brewer(palette = palette) + theme_classic() + ggtitle(paste(GENE, "gene expression")) + 
              theme(axis.title.y = element_text(margin = margin(r = 10), size = title.size), axis.title.x=element_blank()) +
              theme(axis.text = element_text(size = axis.size)) + theme(title = element_text(size = title.size)) + 
              theme(legend.text = element_text(size = legend.text))
            
          } else {
            p <- p + theme_classic() + ggtitle(paste(GENE, "gene expression")) + 
            theme(axis.title.y = element_text(margin = margin(r = 10), size = title.size), axis.title.x=element_blank()) +
            theme(axis.text = element_text(size = axis.size)) + theme(title = element_text(size = title.size)) + 
              theme(legend.text = element_text(size = legend.text))
          }
      }
      if (TYPE == "ECDF"){
          p <- ggplot(gene, aes(genes, colour= label)) + stat_ecdf() + 
                   theme_classic() + ggtitle(paste(GENE, "cumulative gene expression")) + 
                   xlab(paste(GENE, ' expression')) + ylab('Fraction') + 
            theme(axis.title.y = element_text(margin = margin(r = 10), size = title.size), axis.title.x=element_blank()) +
            theme(axis.text = element_text(size = axis.size)) + theme(title = element_text(size = title.size)) + 
              theme(legend.text = element_text(size = legend.text))

      }
      return (p)
  } else {
    if (TYPE == "Den") {
      
      colfunc<-colorRampPalette(c("skyblue3","aquamarine3","yellow","red"))
      
      if (is.null(sgRNA)) {
        grna <- GetAssayData(RDS, assay = "sgRNA")
        grna <- as.matrix(grna)
        grna <- t(grna)
        grna <- as.data.frame(grna)
        da <- colnames(grna)[1]
        target_sgrna_list <-  paste("sgrna_", da, sep = "")
      } else {
        target_sgrna_list = strsplit(sgRNA, ",")[[1]]
        target_sgrna_list = trimws(target_sgrna_list)
        target_sgrna_list <-  paste("sgrna_", target_sgrna_list, sep = "")
        grna <- FetchData(RDS, target_sgrna_list,slot=SLOT)
      }
        grna$grna <- rowSums(grna)
        grna <- subset(grna, grna > 0)
        grna <- grna["grna"]
        da <- FeaturePlot(RDS, target_sgrna_list[1])
        da_1 <- da$data
        num <- nrow(da_1)/3
        da_1 <- da_1[-4]
        data <- merge(da_1, grna, by = 0, all = FALSE)
        if (nrow(data) >= num) {
          p <- ggplot(data, aes(UMAP_1, UMAP_2, z = data[[5]])) + geom_density2d(aes(colour = stat(level))) + 
            scale_color_gradientn(colours = colfunc(4)) + ggtitle(paste("Density of sgRNAs")) + theme_bw()
          p + theme(axis.text = element_text(size = axis.size)) + theme(title = element_text(size = title.size)) + 
            theme(legend.text = element_text(size = legend.text)) 
        } else {
          data <- merge(da_1, grna, by = 0, all = TRUE)
          data[is.na(data)] <- "0"
          rownames(data) <- data[[1]]
          data <- data[-1]
          colnames(data)[4] <- target_sgrna_list[1]
          data[[4]] <- as.numeric(data[[4]])
          da$data <- data
          da + ggtitle("Density of sgRNAs")
        }
        return (p)
    } else {
      message("Please enter a correct type of plot: Dis, Vln, Den, ECDF")
    }
  }
  }
  return (p)
}
