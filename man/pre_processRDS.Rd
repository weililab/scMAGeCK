\name{pre_processRDS}
\alias{pre_processRDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate the imformation of sgRNA into RDS file for the further analysis.
}
\description{
Pre-process the sgRNA count from previous step, and generate the sgRNA expression matrix.
}
\usage{
pre_processRDS(RDS, BARCODE)
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
}
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    Demo <- pre_processRDS(BARCODE = BARCODE, RDS = RDS)
}
