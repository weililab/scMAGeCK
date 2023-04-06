\name{scmageck_best_lambda}
\alias{scmageck_best_lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~
Find the best lambda using negative controls
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function identifies the best value for lambda, a tuning parameter used in a modeling approach that uses negative controls. Negative controls are samples or features that should not have any meaningful association with the outcome of interest. This approach can be useful in situations where there may be confounding factors or batch effects that can affect the results of a model.
}
\usage{
scmageck_best_lambda(
  rds_object,
  bc_frame,
  non_target_ctrl = "NT",
  lambda_seq = 10^seq(-3, 3, length = 100),
  pseudogene_label = "CD47"
  pseudogene_num = 250
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds_object}{
  Seurat object or local RDS file path that contains the scRNA-seq dataset. Alternatively, you can provide the path to an RDS file. The dataset should only contain cells expressing non-targeting gRNA, which will be used as negative controls for lambda tuning.
  }
  \item{bc_frame}{
  A txt file or corresponding data.frame that includes cell identity information generated from the cell identity collection step. The bc_frame should have the same number of rows as the rds_object input, and each row should contain a unique cell barcode or identifier along with its corresponding cell identity.
  }
  \item{non_target_ctrl}{
  A character vector or string containing a comma-separated list of genes that serve as negative controls.
  }
  \item{lambda_seq}{
  A sequence of lambda values to test. Defaults to 10^seq(-3, 3, length = 100).
  }
  \item{pseudogene_label}{
  A character string containing the name of any gene to be used as a pseudogene to evaluate the false positive rate. The pseudogene should be present in the dataset. The defalut value for pseudogene_label is CD47.
  }
  \item{pseudogene_num}{
  An integer value indicating the number of cells that should be selected as pseudogenes for false positive rate evaluation. The default value for pseudogene_num is 250.
  }
}
\value{
The output of the function is a data frame containing the lambda value and the corresponding false positive rate on the negative control samples. Additionally, the function generates a plot showing the relationship between lambda and the false positive rate.
}
\seealso{
\code{\link{scmageck_eff_estimate}} \code{\href{https://bitbucket.org/weililab/scmageck/src/master/}{scMAGeCK BitBucket}}
}
