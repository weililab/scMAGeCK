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
pre_processRDS(BARCODE, RDS, normalize = TRUE, scale = TRUE) 
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
  \item{normalize}{
Whether to perform normalization on sgRNA count matrix 
}
  \item{scale}{
Whether to scale the normalized sgRNA count matrix 
}
}
\value{
Returns an updated Seurat object, with the following modifications:

	- An added "sgrna_guide" assay that contains the normalized and scaled guide counts;
	- An added "sgrna" assay that contains the normalized and scaled target gene counts. Counts from different sgRNAs targeting the same gene are merged;
        - Several columns added to metadata to describe guide assignments (which can be overwritten by calling assign_cell_identity function later).
}
\examples{
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    Demo <- pre_processRDS(BARCODE = BARCODE, RDS = RDS)
}
