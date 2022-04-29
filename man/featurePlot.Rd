\name{featurePlot}
\alias{featurePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect the sgRNA distribution and generate Vlnplot to identity gene regulation between
different cells.
}
\description{
identify how many sgRNAs entered the cells, visualize the counts of sgRNAs distribution and
gene regulation.
}
\usage{
featurePlot(RDS, TYPE = plot.type, BARCODE = NULL, sgRNA = NULL, GENE = NULL, CONTROL = NULL, GROUP2=NULL,SLOT='data',palette = NULL, label.size = 3, axis.size = 12, title.size = 15, legend.text = 10, fill = "#56B4E9")

plot.type
# c("Dis", "Vln", "Den")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDS}{
RDS object from the pre-processRDS step
}
  \item{TYPE}{
Type of the plot.
}
  \item{BARCODE}{
A txt file to include cell identity information, generated from the cell 
identity collection step.
}
 \item{sgRNA}{
Generate whole sgRNAs distribution when sgRNA = NULL, add sgRNAs to see the specific 
sgRNA distribution. Mutiple sgRNAs can be provided, separated by ",". For example, "APC,TP53".
}
 \item{GENE}{
Genes whose expressions are to be compared under different cell coditions. 
Mutiple genes can be provided, as a vector or as a string separated by ",". 
For example, "APC,TP53". 
when provide mutiple genes, it would show the average gene expression. 
}
 \item{CONTROL}{
Set up the sepecifc clusters to compare the gene expression. it would compare gene
expression across the dataset when CONTROL = NULL.
}
 \item{GROUP2}{
Set up the gene name as the control group. it would compare gene
expression across the dataset when GROUP2= NULL.
}
 \item{SLOT}{
Use the slot in Seurat object for plot. May choose 'data', 'scale.data' or 'count'. See GetAssayData funciton in Seurat.
}
  \item{palette}{
The color palette to change the color of VlnPlot.
}
  \item{label.size}{
Text size of label.
}
  \item{axis.size}{
Text size of axis.text.
}
  \item{title.size}{
Text size of axis/pics' title.
}
  \item{legend.text}{
Text size of figure legend.
}
  \item{fill}{
Fill colour.
}
}
\examples{
    
    ### Loading required package
    require(Seurat)
    
    ### BARCODE file contains cell identity information, generated from the cell identity collection step
    BARCODE <- system.file("extdata","barcode_rec.txt",package = "scMAGeCK")
    
    ### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
    RDS <- system.file("extdata","singles_dox_mki67_v3.RDS",package = "scMAGeCK")
    
    ### For using the featurePlot function, it needs to do the preprocessRDS first
    Demo <- pre_processRDS(BARCODE = BARCODE, RDS = RDS)
    
    ### For the sgRNA distribution
    featurePlot(BARCODE = BARCODE, RDS = Demo, TYPE = "Dis")
    
    ### For the density of sgRNA, clustering needed to be done first.
    Demo <- RunUMAP(Demo, dims = 1:10)
    featurePlot(RDS = Demo, sgRNA = NULL, TYPE = "Den")
    
    ###For the Vlnplot to display the gene regulation, take MKI67 for examples
    featurePlot(RDS = Demo, GENE = "MKI67", sgRNA = "TP53", TYPE = "Vln")
}
