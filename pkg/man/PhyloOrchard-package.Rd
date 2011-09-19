\name{PhyloOrchard-package}
\alias{PhyloOrchard-package}
\alias{PhyloOrchard}
\docType{package}
\title{
Important and/or useful phylogenetic datasets
}
\description{
Description: This has trees, source data, and character data of broad utility for phylogenetics. Note that it deliberately does not duplicate much data already in TreeBase. Make sure to cite the original source of the tree/data!
}
\details{
\tabular{ll}{
Package: \tab PhyloOrchard\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-09-19\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
You will typically use data(dataset.name) to load the data. Trees will be in dataset.name$trees. 
}
\author{
Brian C. O'Meara

Maintainer:
omeara.brian@gmail.com
}
\references{
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
}
\examples{
data(BinindaEmondsEtAl2007)
#best tree as default phylo4 format:
library(phylobase)
best.phylo4<-BinindaEmondsEtAl2007$trees$best

#best tree in ape (phylo) format:
library(ape)
best.phylo<-as(best.phylo4,"phylo")

#plot the tree
plot(best.phylo4,show.node.label=TRUE)
}
