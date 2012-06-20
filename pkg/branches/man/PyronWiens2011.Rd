\name{PyronWiens2011}
\alias{PyronWiens2011}
\docType{data}
\title{
Pyron and Wiens 2011 tree of amphibians
}
\description{
This has a congruified timetree (PyronWiens2011) from this paper. Time-scaling used \code{\link{PhyloOrchard}[AlfaroEtAl2009]} as the donor tree and \pkg{treePL} to smooth rates. The tree is in phylo (\pkg{ape}) format.
}
\usage{data(PyronWiens2011)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Pyron, R. A. and J. J. Wiens. 2011. A large-scale phylogeny of Amphibia including over 2800 species, and a revised classification of extant frogs, salamanders, and caecilians. Molecular Phylogenetics and Evolution 61: 543-583.
}

\examples{
#not run
#data(PyronWiens2011)
#amphibianTree<-PyronWiens2011_tree
#plot the tree
#plot(amphibianTree,show.node.label=TRUE)
}
\keyword{datasets}
