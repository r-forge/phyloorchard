\name{HeathEtAl2012}
\alias{HeathEtAl2012}
\docType{data}
\title{
Heath et al. 2012 trees of cute-looking mouse lemurs.
}
\description{
This is a set of 20,000 phylogenetic trees in phylo (ape) format. 
}
\usage{data(HeathEtAl2012)}
\details{
Node age estimates were sampled 
from 6 independent MCMC analyses using DPPDiv (Heath et al. 2012; 
http://cteg.berkeley.edu/software.html#dppdiv) on a fixed topology, using the data and tree from 
Yang and Yoder (2003).
}
\source{
Heath, T. A., M. T. Holder and J. P. Huelsenbeck. 2012. A Dirichlet process prior for estimating 
lineage-specific substitution rates. Mol. Biol. Evol. 29: 939-955.

Yang, Z. and A. D. Yoder. 2003. Comparison of likelihood and Bayesian methods for estimating 
divergence times using multiple gene loci and calibration points, with application to a radiation 
of cute-looking mouse lemur species. Syst. Biol. 52: 705-716. 
}

\examples{
#not run
#data(HeathEtAl2012)

#library(ape)
#treeOne<-HeathEtAl2012$trees[[1]]

#plot the tree
#plot(treeOne,show.node.label=TRUE)

}
\keyword{datasets}
