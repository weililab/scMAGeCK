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
PERMUTATION = NULL, SAVEPATH = "./",LAMBDA=0.01,GENE_FRAC=0.01)
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
that have expressions in at least that fractions of cells are kept. Default: 0.01.
}
}
