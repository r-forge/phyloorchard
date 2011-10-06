\name{BinindaEmondsEtAl2007}
\alias{BinindaEmondsEtAl2007}
\docType{data}
\title{
Bininda-Emonds et al. 2007 mammal supertree (corrected)
}
\description{
This has the best (BinindaEmondsEtAl2007$trees$best), upper (BinindaEmondsEtAl2007$trees$upper), and lower (BinindaEmondsEtAl2007$trees$lower) estimates of branch lengths from this paper. Note that these trees are the corrected ones, not the original ones. The trees are in phylo4 format.
}
\usage{data(BinindaEmondsEtAl2007)}
\details{
Make sure to cite the paper below when you use this dataset.
}
\source{
Olaf R. P. Bininda-Emonds, Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman & Andy Purvis. "The delayed rise of present-day mammals"
Nature 446, 507-512(29 March 2007). doi:10.1038/nature05634}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
#not run
#data(BinindaEmondsEtAl2007)
#best tree as default phylo4 format:
#library(phylobase)
#best.phylo4<-BinindaEmondsEtAl2007$trees$best

#best tree in ape (phylo) format:
#library(ape)
#best.phylo<-as(best.phylo4,"phylo")

#plot the tree
#plot(best.phylo4,show.node.label=TRUE)
}
\keyword{datasets}
