\name{selectPlot}
\alias{selectPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate the selection plot
}
\description{
detect the gene regulation relationship between genes and perturbation by using 
RRA or LR test.
}
\usage{
selectPlot(GENE = NULL, lr_result = NULL, CUTOFF = 0.05, ADJ = "fdr", 
RRA_re1 = NULL, RRA_re2 = NULL, TYPE = select.type, QUALITY = 10)

select.type
#c("lr", "rra")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GENE}{
Genes whose expressions are to be tested under the LR test. Multiple genes can be provided, 
separated by ",". For example, "MKI67,TP53"
}
  \item{lr_result}{
The result from the scmageck-lr step.
}
  \item{CUTOFF}{
Determine the significant pvalue.
}
  \item{ADJ}{
P.adjust.methods. Choose one of correction method.
c("holm","hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
}
  \item{RRA_re1}{
RRA result from the scmageck-rra step.
}
  \item{RRA_re2}{
Optional input. The second RRA result from the scmageck-rra step. Add this input to visualize
gene selection under two different cell condition.
}
  \item{TYPE}{
	  The type of the scMAGeCK results. Can be either "rra" or "lr".
}
  \item{QUALITY}{
The number of single-cells that passes the threshold when use the RRA test, default is 10.
Lower quality could improve the sensitivity but reduce accuracy.
}
}
\examples{
    ### by using RRA test, take MKI67 for example
    ### only works if you have RRA installed
    #selectPlot(RRA_re1 = rra_result, CUTOFF = 0.05, QUALITY = 10, ADJ = "fdr", TYPE = "rra")
    
    ### by using LR test, take MKI67 for example
    
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', SIGNATURE = NULL,
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000, SAVEPATH=NULL, LAMBDA=0.01)
    selectPlot(GENE = "MKI67", lr_result = lr_result, CUTOFF = 0.05, ADJ = "fdr", TYPE = "lr")
}
