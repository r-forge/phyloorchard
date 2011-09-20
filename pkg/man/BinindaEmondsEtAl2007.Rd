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
\format{
  The format is:
List of 1
 $ trees:List of 3
  ..$ best :Formal class 'phylo4' [package "phylobase"] with 6 slots
  .. .. ..@ edge       : int [1:6618, 1:2] 4511 4512 4513 4513 4512 4511 4514 4515 4516 4517 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:2] "ancestor" "descendant"
  .. .. ..@ edge.length: Named num [1:6618] 102.6 53 10.6 10.6 63.6 ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "4511-4512" "4512-4513" "4513-1" "4513-2" ...
  .. .. ..@ label      : Named chr [1:6618] "Tachyglossus_aculeatus" "Zaglossus_bruijni" "Ornithorhynchus_anatinus" "Anomalurus_beecrofti" ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "1" "2" "3" "4" ...
  .. .. ..@ edge.label : Named chr(0) 
  .. .. .. ..- attr(*, "names")= chr(0) 
  .. .. ..@ order      : chr "unknown"
  .. .. ..@ annote     : list()
  ..$ upper:Formal class 'phylo4' [package "phylobase"] with 6 slots
  .. .. ..@ edge       : int [1:6618, 1:2] 4511 4512 4513 4513 4512 4511 4514 4515 4516 4517 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:2] "ancestor" "descendant"
  .. .. ..@ edge.length: Named num [1:6618] 114 43.6 8.6 8.6 52.2 24.4 48.6 0.1 0.2 2.9 ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "4511-4512" "4512-4513" "4513-1" "4513-2" ...
  .. .. ..@ label      : Named chr [1:6618] "Tachyglossus_aculeatus" "Zaglossus_bruijni" "Ornithorhynchus_anatinus" "Anomalurus_beecrofti" ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "1" "2" "3" "4" ...
  .. .. ..@ edge.label : Named chr(0) 
  .. .. .. ..- attr(*, "names")= chr(0) 
  .. .. ..@ order      : chr "unknown"
  .. .. ..@ annote     : list()
  ..$ lower:Formal class 'phylo4' [package "phylobase"] with 6 slots
  .. .. ..@ edge       : int [1:6618, 1:2] 4511 4512 4513 4513 4512 4511 4514 4515 4516 4517 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:2] "ancestor" "descendant"
  .. .. ..@ edge.length: Named num [1:6618] 91.2 62.5 12.5 12.5 75 13.1 45 8.3 1.4 4.6 ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "4511-4512" "4512-4513" "4513-1" "4513-2" ...
  .. .. ..@ label      : Named chr [1:6618] "Tachyglossus_aculeatus" "Zaglossus_bruijni" "Ornithorhynchus_anatinus" "Anomalurus_beecrofti" ...
  .. .. .. ..- attr(*, "names")= chr [1:6618] "1" "2" "3" "4" ...
  .. .. ..@ edge.label : Named chr(0) 
  .. .. .. ..- attr(*, "names")= chr(0) 
  .. .. ..@ order      : chr "unknown"
  .. .. ..@ annote     : list()
}
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
