ncbiTaxonomy<-function(names.desired=c("scientific","common"),mimimum.rank=c("any","species","genus","family","order","class","phylum","kingdom"))) {
  library("R.utils")
  system("mkdir ncbi_taxdmp")
  setwd("ncbi_taxdmp")
  pathname <- downloadFile("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
  system(paste(" tar -xzf ",pathname))
  names.dmp<-read.table("names.dmp",header=FALSE, sep="|",strip.white=TRUE,fill=TRUE) #if there is a "#" everything after that isn't read
  names.dmp<-names.dmp[,1:4]
  names(names.dmp)<-c("tax_id", "name_txt", "unique_name", "name_class")
  nodes.dmp<-read.table("nodes.dmp",header=FALSE, sep="|",strip.white=TRUE,fill=TRUE)
  nodes.dmp<-nodes.dmp[,c(1:5,11:13)]
  names(nodes.dmp)<-c("tax_id","parent_tax_id","rank","embl_code","division_id","GenBank_hidden_flag","hidden_subtree_root_flag","comments")
  
  names.pruned<-c()
  if(names.desired=="scientific") {
     names.pruned<-subset(names.dmp,names.dmp[,4]=="scientific name")
  }
  else if(names.desired=="common") {
     names.pruned<-subset(names.dmp,names.dmp[,4]=="genbank common name")
  }
  
  tree.matrix.raw<-merge(nodes.dmp,names.pruned,all.x=TRUE)
  
  system("rm -r ncbi_taxdmp")
}