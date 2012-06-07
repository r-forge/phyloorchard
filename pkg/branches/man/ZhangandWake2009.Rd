\name{ZhangandWake2009}
\alias{ZhangandWake2009}
\docType{data}
\title{
Zhang and Wake 2009 tree of caecilians. This tree is a bit of an outlier here as it is not that big and also not species-level.
}
\description{
This has a single timetree (ZhangandWake2009$trees$best) from this paper. The tree is in phylo4 format.
}
\usage{data(ZhangandWake2009)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Zhang, P. and M. H. Wake. 2009. Mitogenomic perspective on the phylogeny and biogeography of living caecilians (Amphibia:Gymnophiona).  Mol.  Phyl. Evol. 53:479-491. 
}

\examples{
#not run
#data(ZhangandWake2009)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-ZhangandWake2009$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
