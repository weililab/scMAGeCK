\name{selectPlot}
\alias{selectPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate the selection plot
}
\description{
detect the gene regulation relationship between genes and perturbation by using 
RRA or LR test.
}
\usage{
selectPlot(GENE = NULL, lr_result = NULL, CUTOFF = 0.05, ADJ = "fdr", 
RRA_re1 = NULL, RRA_re2 = NULL, TYPE = select.type, QUALITY = 10)

select.type
#c("lr", "rra")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GENE}{
Genes whose expressions are to be tested under the LR test. Multiple genes can be provided, 
separated by ",". For example, "MKI67,TP53"
}
  \item{lr_result}{
The result from the scmageck-lr step.
}
  \item{CUTOFF}{
Determine the significant pvalue.
}
  \item{ADJ}{
P.adjust.methods. Choose one of correction method.
c("holm","hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
}
  \item{RRA_re1}{
RRA result from the scmageck-rra step.
}
  \item{RRA_re2}{
Optional input. The second RRA result from the scmageck-rra step. Add this input to visualize
gene selection under two different cell condition.
}
  \item{QUALITY}{
The number of single-cells that passes the threshold when use the RRA test, default is 10.
Lower quality could improve the sensitivity but reduce accuracy.
}
}
\examples{
    ### by using RRA test, take MKI67 for example
    selectPlot(RRA_re1 = rra_result, CUTOFF = 0.05, QUALITY = 10, ADJ = "fdr", TYPE = "rra")
    
    ### by using LR test, take MKI67 for example
    selectPlot(GENE = "MKI67", lr_result = lr_result, CUTOFF = 0.05, ADJ = "fdr", TYPE = "lr")
}