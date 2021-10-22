\name{scmageck_lr}
\alias{scmageck_lr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use linear regression to test the association of gene knockout with all possible genes
}
\description{
echo "Use linear regression to test the association of gene knockout with all possible genes"
}
\usage{
scmageck_lr(BARCODE, RDS, NEGCTRL, SELECT_GENE=NULL, LABEL = NULL,  
PERMUTATION = NULL, SIGNATURE = NULL, SAVEPATH = "./",LAMBDA=0.01,GENE_FRAC=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BARCODE}{
A txt file to include cell identity information, generated from the cell 
identity collection step.
}
  \item{RDS}{
A Seurat object or local RDS file path that contains the scRNA-seq dataset. 
Note that the dataset has to be normalized and scaled.
}
  \item{NEGCTRL}{
The name of the genes (separated by ",") served as negative controls.
}
  \item{SELECT_GENE}{
The list of genes for regression. By default, all genes in the table are 
subject to regression. 
}
  \item{LABEL}{
The label of the output file.
}
  \item{PERMUTATION}{
The number of permutations for p value calculation.
}
  \item{SAVEPATH}{
The save path of result. Default save path is the current working directory. 
If you don't need save the result, set SAVEPATH as NULL.
}
  \item{LAMBDA}{
A paramter for the LR model for ridge regression. Default: 0.01.
}
  \item{GENE_FRAC}{
A paramter for filtering low expressed genes. By default, only genes 
that have expressions in at least that fractions of cells (in raw count table) are kept. 
If raw count table is not found, will look into scaled expression, 
and only keep genes whose expression is greater than zero in at least that fraction of cells.
Default: 0.01.
}
  \item{SIGNATURE}{
A GMT text file, the format must be as follows:(1)Column 1: name of gene set; 
(2)Colum 2: Empty, or the information about gene set e.g. the source of the gene set; 
(3)Column 3 and onwards: ids of genes beloging to a particular gene set. 
Note that this parameter for applying LR model for the gene set analysis.
Refercence:http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
}
}
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', SIGNATURE = NULL,
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000, SAVEPATH=NULL, LAMBDA=0.01)
    lr_score <- lr_result[1]
    lr_score_pval <- lr_result[2]
    head(lr_score_pval)
}
\value{
The result for object \code{RDS}
}
