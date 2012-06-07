\name{BergmannEtAl2012}
\alias{BergmannEtAl2012}
\docType{data}
\title{
Bergmann and Irschick 2012 tree of squamates
}
\description{
This has a single timetree (BergmannEtAl2012$trees$best) from this paper. The tree is in phylo4 format.
}
\usage{data(BergmannEtAl2012)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Bergmann, P. J. and D. J. Irschick. Vertebral evolution and the diversification of squamate reptiles. Evolution 66: 1044-1058.
}

\examples{
#not run
#data(BergmannEtAl2012)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-BergmannEtAl2012$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
