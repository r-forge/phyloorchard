\name{EastmanEtAlUnpublished}
\alias{EastmanEtAlUnpublished}
\docType{data}
\title{
Eastman et al. unpublished tree of salamanders
}
\description{
This has a single timetree of salamanders (EastmanEtAlUnpublished$trees$best). The tree is in phylo4 format.
}
\usage{data(EastmanEtAlUnpublished)}
\details{
Tree is unpublished, contact Jon Eastman for more information}
\source{

}

\examples{
#not run
#data(EastmanEtAlUnpublished)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-EastmanEtAlUnpublished$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
