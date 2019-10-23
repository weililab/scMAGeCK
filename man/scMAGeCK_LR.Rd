\name{scMAGeCK_LR}
\alias{scMAGeCK_LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use linear regression to test the association of gene knockout with all possible genes
}
\description{
echo "Use linear regression to test the association of gene knockout with all possible genes"
}
\usage{
scMAGeCK_LR(BARCODE, RDS, NEGCTRL, SELECT_GENE=NULL, LABEL = NULL, PERMUTATION = NULL, SAVEPATH = "./",LAMBDA=0.01,GENE_FRAC=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BARCODE}{
A txt file to include cell identity information, generated from the cell identity collection step.
}
  \item{RDS}{
A Seurat object or local RDS file path that contains the scRNA-seq dataset. Note that the dataset has to be normalized and scaled.
}
  \item{NEGCTRL}{
The name of the genes (separated by ",") served as negative controls.
}
  \item{SELECT_GENE}{
The list of genes for regression. By default, all genes in the table are subject to regression. 
}
  \item{LABEL}{
The label of the output file.
}
  \item{PERMUTATION}{
The number of permutations for p value calculation.
}
  \item{SAVEPATH}{
The save path of result. Default save path is the current working directory. If you don't need save the result, set SAVEPATH as NULL.
}
  \item{LAMBDA}{
A paramter for the LR model.
}
  \item{GENE_FRAC}{
A paramter for filtering low expressed genes.
}
}
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    lr_result <- scMAGeCK_LR(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', 
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000, SAVEPATH=NULL, LAMBDA=0.01)
    lr_score <- lr_result[1]
    lr_score_pval <- lr_result[2]
    head(lr_score_pval)
}
\value{
The result for object \code{RDS}
}
