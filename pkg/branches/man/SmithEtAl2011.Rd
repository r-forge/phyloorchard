\name{SmithEtAl2011}
\alias{SmithEtAl2011}
\docType{data}
\title{
Smith et al. 2011 angiosperm tree
}
\description{
This has the best (SmithEtAl2011$trees$best) likelihood topology. The tree is in phylo4 format.
}
\usage{data(SmithEtAl2011)}
\details{
Make sure to cite the paper and dryad reference below when you use this dataset.
}
\source{
Smith SA, Beaulieu JM, Stamatakis A, Donoghue MJ (2011) Understanding angiosperm diversification using small and large phylogenetic trees. American Journal of Botany, online in advance of print. doi:10.3732/ajb.1000481

Smith SA, Beaulieu JM, Stamatakis A, Donoghue MJ (2011) Data from: Understanding angiosperm diversification using small and large phylogenetic trees. Dryad Digital Repository. doi:10.5061/dryad.8790
\examples{
#not run
#data(SmithEtAl2011)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-SmithEtAl2011$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
