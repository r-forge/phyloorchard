\name{TreeBase}
\alias{treebase}
\alias{Treebase}
\alias{TreeBASE}
\docType{data}
\title{
TreeBase
}
\description{
This has the the trees from TreeBase as of Oct. 5, 2011. They were downloaded using Carl Boettiger's rtreebase package.
}
\usage{data(TreeBase)}
\details{
This has the tree files from TreeBase. Note that they are stored in phylo (\pkg{ape}) format.  The files are stored in a list. Each element in the list is one study. Each element of the list is itself a list, containing information such as the study id, year of publication, authors, and other citation information. You should cite the authors of the study, especially if a single tree is being reused (as opposed to, say, looking at the distribution of tree balance across thousands of studies). The trees themselves are in the $trees object in each list: TreeBase[[1]]$trees has the MultiPhylo object containing tree(s) from the first study. 
}
\source{
TreeBase: http://www.treebase.org/
}
\references{
}
\examples{
#not run
#data(TreeBase)
#plot the first tree
#plot(TreeBase[[1]]$trees[[1]])

#put the first tree in phylo4 format
#tree1<-as(TreeBase[[1]]$trees[[1]],"phylo4")

}
\keyword{datasets}
