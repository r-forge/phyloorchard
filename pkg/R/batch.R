source("ncbiTaxonomy.R")
result<-ncbiTaxonomy(names.desired="scientific",minimum.rank="species")
save(result,file="NCBIresult.RData",compress=TRUE)
