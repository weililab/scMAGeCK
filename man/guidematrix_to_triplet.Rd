\name{guidematrix_to_triplet}
\alias{guidematrix_to_triplet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert guide matrix into a dataframe containing cell, barcode, read and UMI count.
}
\description{
Convert guide matrix into a dataframe, where columns include cell, barcode, 
read count and UMI count.
}
\usage{
guidematrix_to_triplet(count_mat, RDS) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count_mat}{
A guide count matrix where columns are cell names and rows are guide names
} 
  \item{RDS}{
A Seurat object or local RDS file path that contains the scRNA-seq dataset. 
Note that the dataset has to be normalized and scaled.
}
}
\value{
  A data frame that contains cell names, guide names, read count and UMI count as column.
}
\examples{
}
