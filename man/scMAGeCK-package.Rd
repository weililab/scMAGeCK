\name{scMAGeCK-package}
\alias{scMAGeCK-package}
\alias{scmageck}
\docType{package}
\title{
\packageTitle{scMAGeCK}
}
\description{
\packageDescription{scMAGeCK}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{scMAGeCK}
\packageIndices{scMAGeCK}
    scMAGeCK is a computational model to identify genes associated with 
    multiple expression phenotypes from CRISPR screening coupled with 
    single-cell RNA sequencing data (CROP-seq).scMAGeCK is based on our 
    previous MAGeCK and MAGeCK-VISPR models for pooled CRISPR screens. 
    
    The scMAGeCK manuscript can be found at bioRxiv(https://www.biorxiv.org/content/10.1101/658146v1/).
}
\author{
\packageAuthor{scMAGeCK}

Maintainer: \packageMaintainer{scMAGeCK}
}

\keyword{ package }
\examples{
    ### BARCODE file contains cell identity information, generated from 
    ### the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    ### RDS can be a Seurat object or local RDS file path that contains 
    ### the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    ### Set RRA executable file path. 
    ###   You can generate RRA executable file by following commands:
    ###   wget https://bitbucket.org/weililab/scmageck/downloads/RRA_0.5.9.zip
    ###   unzip RRA_0.5.9.zip
    ###   cd RRA_0.5.9
    ###   make
    RRAPATH <- "/Library/RRA_0.5.9/bin/RRA"
    
    target_gene <- "MKI67"
    
    rra_result <- scmageck_rra(BARCODE=BARCODE, RDS=RDS, GENE=target_gene, 
                               RRAPATH=RRAPATH, LABEL='dox_mki67', 
                               NEGCTRL=NULL, KEEPTMP=FALSE, 
                               PATHWAY=FALSE, SAVEPATH=NULL)
    head(rra_result)
    
    lr_result <- scmageck_lr(BARCODE=BARCODE, RDS=RDS, LABEL='dox_scmageck_lr', 
            NEGCTRL = 'NonTargetingControlGuideForHuman', PERMUTATION = 1000,
            SAVEPATH=NULL, LAMBDA=0.01)
    lr_score <- lr_result[1]
    lr_score_pval <- lr_result[2]
    head(lr_score_pval)
}
