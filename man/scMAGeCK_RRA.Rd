\name{scmageck_rra}
\alias{scmageck_rra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use RRA to test the association of gene knockout with certain marker expression
}
\description{
echo "Use RRA to test the association of gene knockout with certain marker expression"
}
\usage{
scmageck_rra(BARCODE, RDS, GENE, RRAPATH = NULL, LABEL = NULL, NEGCTRL = NULL, 
KEEPTMP = FALSE, PATHWAY = FALSE, SAVEPATH = "./")
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
  \item{GENE}{
Genes whose expressions are to be tested. Multiple genes can be provided, 
separated by ",". For example, "MKI67,TP53"
}
  \item{RRAPATH}{
The path to the RRA program, if RRA cannot be found in the PATH environment variable.
}
  \item{LABEL}{
The label of the output file.
}
  \item{NEGCTRL}{
The name of the negative control gene. For example, "NonTargetingControlGuideForHuman". 
Default is NULL (do not use any negative controls).
}
  \item{KEEPTMP}{
Keep temporary files.
}
  \item{PATHWAY}{
Treat genes in --GENE option as a pathway. In other words, the averaged 
expression of these genes will be used for testing.
}
  \item{SAVEPATH}{
The save path of result. Default save path is the current working directory. 
If you don't need save the result, set SAVEPATH as NULL.
}
}
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    target_gene <- "MKI67"
    
    rra_result <- scmageck_rra(BARCODE=BARCODE, RDS=RDS, GENE=target_gene,
             LABEL='dox_mki67', NEGCTRL=NULL, KEEPTMP=FALSE, PATHWAY=FALSE, SAVEPATH=NULL)
    head(rra_result)
}
\value{
The result for object \code{RDS}
}

