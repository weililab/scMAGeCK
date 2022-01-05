read_gmt_file <- function(gmt_path) {
   x <- scan(gmt_path, what = "", sep = "\n")
    x <- strsplit(x, "\t") # split string by white space
    max.col <- max(sapply(x, length))
    cn <- paste("V", 1:max.col, sep = "")
    gmt <- read.table(gmt_path, fill = TRUE, col.names = cn)
    # gmt <- read.delim(gmt_path, header = FALSE)
    gmt <- t(as.matrix(gmt))
    colnames(gmt) <- gmt[1, ]
    gmt <- gmt[-1:-2, ]
    message(paste("Total signature records:", ncol(gmt)))
    return (gmt)
}
