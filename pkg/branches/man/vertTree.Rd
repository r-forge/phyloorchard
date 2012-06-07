\name{vertTree}
\alias{vertTree}
\docType{data}
\title{
Giant glomogram of vertebrates. Use at your own risk.
}
\description{
This has a single timetree of ~15,000 vertebrates (vertTree$trees$best). The tree is in phylo4 format.
}
\usage{data(vertTree)}
\details{
Tree is unpublished, contact Jon Eastman for more information}
\source{

}

\examples{
#not run
#data(vertTree)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-vertTree$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
