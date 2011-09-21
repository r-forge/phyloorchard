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
  
  tree.matrix.ncbi<-merge(nodes.dmp,names.pruned,all.x=TRUE)
  tree.matrix.ncbi[,1]<-as.integer(tree.matrix.ncbi[,1])
  tree.matrix.ncbi[,2]<-as.integer(tree.matrix.ncbi[,2])

  ancestors.ncbi<-as.integer(levels(as.factor(tree.matrix.ncbi$parent_tax_id)))
  descendants.all.ncbi<-as.integer(levels(as.factor(tree.matrix.ncbi$tax_id)))
  tips.ncbi<-descendants.all.ncbi[is.na(match(descendants.all.ncbi,ancestors.ncbi))]
  if (minimum.rank!="any") {
    tips.ncbi<-as.integer(levels(as.factor(tree.matrix.ncbi[which(tree.matrix.ncbi$rank==minimum.rank),]$tax_id)))
  }
  tips.translation<-data.frame(tips.phylo=c(1:length(tips.ncbi)),tips.ncbi=tips.ncbi)
  tips.name<-as.character(tree.matrix.ncbi$name_txt[match(tips.ncbi,tree.matrix.ncbi$tax_id)])
  tips.name[is.na(tips.name)]<-""
  all.nodes.translation<-tips.translation
  names(all.nodes.translation)<-c("nodes.phylo","nodes.ncbi")
  phylo.matrix<-matrix(,nrow=1,ncol=2)
  for (tip.phylo in 1:length(tips.ncbi)) {
    continue=TRUE
    current.node.phylo<-tip.phylo
    current.node.ncbi<-tips.translation[which(tips.translation$tips.phylo==current.node.phylo),]$tips.ncbi
    
    while(continue) {
      print(paste("current.node.phylo=",current.node.phylo))
      print(paste("current.node.ncbi=",current.node.ncbi))

      parent_tax_id.ncbi<-tree.matrix.ncbi[which(tree.matrix.ncbi$tax_id==current.node.ncbi),]$parent_tax_id
      print(paste("parent_tax_id.ncbi=",parent_tax_id.ncbi))

      if(is.na(match(parent_tax_id.ncbi,all.nodes.translation$nodes.ncbi))) {
        print(paste("parent_tax_id.ncbi=",parent_tax_id.ncbi))
        print(c(1+max(c(length(tips.ncbi),max(all.nodes.translation$nodes.phylo,na.rm=TRUE)),na.rm=TRUE),parent_tax_id.ncbi))
        new.df<-data.frame(nodes.phylo=1+max(c(length(tips.ncbi),max(all.nodes.translation$nodes.phylo,na.rm=TRUE)),na.rm=TRUE),nodes.ncbi=parent_tax_id.ncbi)
        print(new.df)
        dput(new.df)
        dput(all.nodes.translation)
        all.nodes.translation<-rbind(all.nodes.translation,new.df)
      }
      print(all.nodes.translation)
      parent_tax_id.phylo<-all.nodes.translation[which(all.nodes.translation$nodes.ncbi==parent_tax_id.ncbi),]$nodes.phylo
      phylo.matrix<-rbind(phylo.matrix,c(parent_tax_id.phylo,current.node.phylo))
      print(c(parent_tax_id.phylo,current.node.phylo))
      current.node.ncbi<-parent_tax_id.ncbi
      current.node.phylo<-parent_tax_id.phylo
      if (!is.na(match(current.node.phylo,phylo.matrix[,2]))) {
        continue=FALSE
      }
    }
  }
  phylo.matrix<-phylo.matrix[-1,] #remove intial values
  phylo.matrix<-phylo.matrix[order(phylo.matrix[,2]),]
  original.root<-phylo.matrix[which(phylo.matrix[,1]==phylo.matrix[,2]),1]
  phylo.matrix[which(phylo.matrix[,1]==phylo.matrix[,2]),1]<-0 #set root to zero for phylobase
  all.nodes.translation[which(all.nodes.translation$nodes.phylo==original.root),1]<-0
  ncbi.phylo4<-phylo4(x=phylo.matrix)
  tipLabels(ncbi.phylo4)<-tips.name
  #system("rm -r ncbi_taxdmp")
}