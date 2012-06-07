\name{Oaks2011}
\alias{Oaks2011}
\docType{data}
\title{
Oakes 2011 tree of crocodylians
}
\description{
This has a single timetree (Oaks2011$trees$best) from this paper. The tree is in phylo4 format.
}
\usage{data(Oaks2011)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Oaks, J. R. 2011. A time-calibrated species tree of crocodylia reveals a recent
radiation of the true crocodyles. Evolution 11: 3285-3297. 
}

\examples{
#not run
#data(Oaks2011)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-Oaks2011$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
