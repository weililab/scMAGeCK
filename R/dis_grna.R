dis_grna <- function(BARCODE, label.size = 3, axis.size = 12, title.size = 15, fill = "#56B4E9"){
  data <- read.delim(BARCODE)
  L <- unique(data$cell)
  output <- NULL
  for (i in L) {
    a <- nrow(subset(data, cell == i))
    m <- matrix(c(i, a), ncol = 2, dimnames = list(i, c("cell", "number_gRNA")))
    output <- rbind(output, m)
  }
  
  output <- as.data.frame(output)
  output$number_gRNA <- as.numeric(output$number_gRNA)
  
  N <- unique(output$number_gRNA)
  output_1 <- NULL
  for (t in N) {
    b <- nrow(subset(output, number_gRNA == t))
    c <- matrix(c(t, b), ncol = 2, dimnames = list(t, c("number_gRNA", "ncell")))
    output_1 <- rbind(output_1, c)
  }
  
  # cell0 <- matrix(c("0", est_cell - count(output)$n), ncol = 2, dimnames = list(0, c("number_gRNA", "ncell")))
  # output_1 <- rbind(cell0, output_1)
  
  output_1 <- as.data.frame(output_1)
  ggplot(output_1, aes(number_gRNA, ncell)) +
    geom_col(width = 0.4, fill = fill, colour = "black") + 
    ggtitle("sgRNA distribution") + geom_text(aes(label = ncell), size = label.size, hjust = 0.5, vjust = .01) + theme_bw() +
    theme(axis.text = element_text(angle = 90, size = axis.size)) + theme(plot.title = element_text(size = title.size)) +
    theme(axis.title = element_text(size = title.size))
}
TRUE