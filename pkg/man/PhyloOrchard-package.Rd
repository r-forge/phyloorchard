\name{PhyloOrchard-package}
\alias{PhyloOrchard-package}
\alias{PhyloOrchard}
\docType{package}
\title{
Chronograms of life
}
\description{
Description: This has chronograms (trees with branch lengths proportional to time). It may have a single tree from a paper, and thus it is in ape's phylo format, or multiple trees, in which case they would be in ape's multiphylo format. In the case of multiple trees, they could be one with median, lower, and upper bounds on branch lengths, or they could be a set of post-burnin BEAST chronograms. Make sure to cite the original source of the tree/data!
}
\details{
\tabular{ll}{
Package: \tab PhyloOrchard\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2012-01-30\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
You will typically use data(dataset.name) to load the data. 
}
\author{
Brian C. O'Meara, Luke Harmon, Jon Eastman

Maintainer:
Brian C. O'Meara <omeara.brian@gmail.com>
}
\keyword{ package }
\seealso{
ape
}
\examples{
#not run
#data(BinindaEmondsEtAl2007)
#library(ape)
#plot(BinindaEmondsEtAl2007[[1]],show.node.label=TRUE)
}
