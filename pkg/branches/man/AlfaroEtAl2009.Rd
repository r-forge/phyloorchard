\name{AlfaroEtAl2009}
\alias{AlfaroEtAl2009}
\docType{data}
\title{
Alfaro et al. 2009 tree of gnathostomes
}
\description{
This has a single timetree (AlfaroEtAl2009$trees$best) from this paper. The tree is in phylo4 format.
}
\usage{data(AlfaroEtAl2009)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Alfaro, M. E., F. Santini, C. Brock, H. Alamillo, A. Dornburg, D. L. Rabosky, G. Carnevale, 
and L. J. Harmon. 2009. Nine exceptional radiations plus high turnover explain species diversity 
in jawed vertebrates. PNAS 106:13410-13414.
}

\examples{
#not run
#data(AlfaroEtAl2009)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-AlfaroEtAl2009$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
