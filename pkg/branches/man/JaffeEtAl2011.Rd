\name{JaffeEtAl2011}
\alias{JaffeEtAl2011}
\docType{data}
\title{
Jaffe et al. 2011 tree of turtles
}
\description{
This has a single timetree (JaffeEtAl2011$trees$best) from this paper. The tree is in phylo4 format.
}
\usage{data(JaffeEtAl2011)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Jaffe, A. L., G. J. Slater, and M. E. Alfaro. The evolution of island gigantism and body size variation in tortoises and turtles. Biology Letters 7: 558-561.
}

\examples{
#not run
#data(JaffeEtAl2011)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-JaffeEtAl2011$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
