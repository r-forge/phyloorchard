\name{AlfaroEtAlUnpublished}
\alias{AlfaroEtAlUnpublished}
\docType{data}
\title{
Alfaro et al. unpublished tree of fish
}
\description{
This has a single timetree of fish (AlfaroEtAlUnpublished$trees$best). The tree is in phylo4 format.
}
\usage{data(AlfaroEtAl2009)}
\details{
Tree is unpublished, contact Mike Alfaro for more information}
\source{

}

\examples{
#not run
#data(AlfaroEtAlUnpublished)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-AlfaroEtAlUnpublished$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
