returnNodesWithJustOneDescendant<-function(phylo.matrix) { 
   singleDescendantNodes<-c()
   factors<-as.factor(phylo.matrix[,1])
   for(i in 1:nlevels(factors)) {
      if (length(which(factors==levels(factors)[i]))==1) {
        singleDescendantNodes<-append(singleDescendantNodes, levels(factors)[i])
      }
   }
   return(singleDescendantNodes)
}


pruneNodesWithOneDescendant<-function(phylo.matrix) {
  #"ancestor first col descendent second col"
  singleDescendantNodes<-returnNodesWithJustOneDescendant(phylo.matrix)
  print(paste("there are ",length(singleDescendantNodes),"nodes with one descendant and ",dim(phylo.matrix)[1],"rows in phylo.matrix"))
  while(length(singleDescendantNodes)>0) {
    badNode<-as.integer(singleDescendantNodes[1])
    ancestorNode<-phylo.matrix[which(phylo.matrix[,2]==badNode),1] #find the ancestor
    if (length(ancestorNode)>0) { #if the root has one descendant, this will be false.
      phylo.matrix[which(phylo.matrix[,1]==badNode),1]<-ancestorNode[1] #point the descendant node to the ancestor
      phylo.matrix<-phylo.matrix[-which(phylo.matrix[,2]==badNode),] #and delete the edge pointing the bad node
    }
    else {
      phylo.matrix<-phylo.matrix[-which(phylo.matrix[,1]==badNode),] #delete the edge descended from the root node
    }
    phylo.matrix[which(phylo.matrix[,1]>badNode),1]<-phylo.matrix[which(phylo.matrix[,1]>badNode),1]-1 #everything has to drop down one in value
    phylo.matrix[which(phylo.matrix[,2]>badNode),2]<-phylo.matrix[which(phylo.matrix[,2]>badNode),2]-1 #everything has to drop down one in value
    singleDescendantNodes<-returnNodesWithJustOneDescendant(phylo.matrix)
    print(paste("there are ",length(singleDescendantNodes),"nodes with one descendant and ",dim(phylo.matrix)[1],"rows in phylo.matrix"))
  }
  return(phylo.matrix)
}

ncbiTaxonomy<-function(names.desired=c("scientific","common"),minimum.rank=c("any","species","genus","family","order","class","phylum","kingdom"),download=TRUE) {
   library("R.utils")
  library("phylobase")
  if(download) {
    system("mkdir ncbi_taxdmp")
    setwd("ncbi_taxdmp")
    print("Downloading from NCBI...")
    pathname <- downloadFile("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    print("Processing downloaded files...")
    system(paste(" tar -xzf ",pathname))
    system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' names.dmp")
   system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' nodes.dmp")
  }
  print("Loading into R...")
  names.dmp<-read.table("names.dmp",header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE) 
  names.dmp<-names.dmp[,1:4]
  names(names.dmp)<-c("tax_id", "name_txt", "unique_name", "name_class")
  nodes.dmp<-read.table("nodes.dmp",header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE)
  nodes.dmp<-nodes.dmp[,c(1:5,11:13)]
  names(nodes.dmp)<-c("tax_id","parent_tax_id","rank","embl_code","division_id","GenBank_hidden_flag","hidden_subtree_root_flag","comments")
  
  names.pruned<-c()
  if(names.desired=="scientific") {
     names.pruned<-subset(names.dmp,names.dmp[,4]=="scientific name")
  }
  if(names.desired=="common") {
     names.pruned<-subset(names.dmp,names.dmp[,4]=="genbank common name")
  }
  print("Files successfully loaded in R")
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
  tips.name[is.na(tips.name)]<-tips.ncbi[is.na(tips.name)]
  all.nodes.translation<-tips.translation
  names(all.nodes.translation)<-c("nodes.phylo","nodes.ncbi")
  phylo.matrix<-matrix(,nrow=1,ncol=5)
  names(phylo.matrix)<-c("ancestor","current","current.ncbi","current.sci.label","current.common.label")
  for (tip.phylo in 1:length(tips.ncbi)) {
    continue=TRUE
    current.node.phylo<-tip.phylo
    current.node.ncbi<-tips.translation[which(tips.translation$tips.phylo==current.node.phylo),]$tips.ncbi
    print(paste("Now doing tip ",tip.phylo,"of ",length(tips.ncbi)))
    #print(paste("current.node.phylo=",current.node.phylo))
    #print(paste("current.node.ncbi=",current.node.ncbi))
    #print ("working")
    while(continue) {
      #cat(".")
      #print(paste("current.node.phylo=",current.node.phylo))
      #print(paste("current.node.ncbi=",current.node.ncbi))
      parent_tax_id.ncbi<-tree.matrix.ncbi[which(tree.matrix.ncbi$tax_id==current.node.ncbi),]$parent_tax_id
      #print(paste("parent_tax_id.ncbi=",parent_tax_id.ncbi))

      if(is.na(match(parent_tax_id.ncbi,all.nodes.translation$nodes.ncbi))) {
        #print(paste("parent_tax_id.ncbi=",parent_tax_id.ncbi))
        print(c(1+max(c(length(tips.ncbi),max(all.nodes.translation$nodes.phylo,na.rm=TRUE)),na.rm=TRUE),parent_tax_id.ncbi))
        new.df<-data.frame(nodes.phylo=1+max(c(length(tips.ncbi),max(all.nodes.translation$nodes.phylo,na.rm=TRUE)),na.rm=TRUE),nodes.ncbi=parent_tax_id.ncbi)
        #print(new.df)
        #dput(new.df)
        #dput(all.nodes.translation)
        all.nodes.translation<-rbind(all.nodes.translation,new.df)
      }
      #print(all.nodes.translation)
      current.node.sci.label<-names.dmp[which(names.dmp$tax_id==current.node.ncbi & names.dmp$name_class=="scientific name"),]$name_txt[1]
      current.node.common.label<-names.dmp[which(names.dmp$tax_id==current.node.ncbi & names.dmp$name_class=="genbank common name"),]$name_txt[1]
      parent_tax_id.phylo<-all.nodes.translation[which(all.nodes.translation$nodes.ncbi==parent_tax_id.ncbi),]$nodes.phylo
      phylo.matrix<-rbind(phylo.matrix,c(parent_tax_id.phylo,current.node.phylo,current.node.ncbi,current.node.sci.label,current.node.common.label))
      print(c(parent_tax_id.phylo,current.node.phylo))
      current.node.ncbi<-parent_tax_id.ncbi
      current.node.phylo<-parent_tax_id.phylo
      if (!is.na(match(current.node.phylo,phylo.matrix[,2]))) {
        continue=FALSE
      }
    }
  }
  phylo.matrix<-phylo.matrix[-1,] #remove intial values
  phylo.matrix<-data.frame(phylo.matrix,stringsAsFactors=FALSE)
  names(phylo.matrix)<-c("ancestor","current","current.ncbi","current.sci.label","current.common.label")

  for (i in 1:3) {
    phylo.matrix[,i]<-as.numeric(phylo.matrix[,i])
  }
  phylo.matrix<-phylo.matrix[order(phylo.matrix[,2]),]
  phylo.matrix<-phylo.matrix[-which(phylo.matrix[,1]==phylo.matrix[,2]),] #get rid of node with itself as an ancestor
  
  phylo.matrix<-pruneNodesWithOneDescendant(phylo.matrix)
  #phylo.matrix.for.ape<-phylo.matrix


  
  original.root.rows<-which(!phylo.matrix[,1] %in% phylo.matrix[,2])
  original.root<-(phylo.matrix[original.root.rows,1])[1]
  tip.rows<-which(!phylo.matrix[,2] %in% phylo.matrix[,1])
  phylo.matrix<-rbind(phylo.matrix,data.frame(ancestor=0,current=original.root,current.ncbi=0,current.sci.label="Life",current.common.label="Life")) #set root to zero for phylobase
  phylo.matrix$current.common.label[which(is.na(phylo.matrix$current.common.label))]<-" "
  phylo.matrix$current.sci.label[which(is.na(phylo.matrix$current.sci.label))]<-" "

  x<-matrix(as.numeric(c(phylo.matrix[,1],phylo.matrix[,2])),ncol=2,byrow=FALSE)
   ncbi.phylo4<-phylo4(x=x)
   labels(ncbi.phylo4,type="tip")<-phylo.matrix$current.sci.label[tip.rows]
   #ncbi.phylo<-reorder(as.phylo(as(ncbi.phylo4,"phylog")))  #has no internal labels
  #ncbi.phylo4.all.labels<-ncbi.phylo4
  #try(labels(ncbi.phylo4.all.labels,type="all")<-phylo.matrix$current.sci.label)
  ncbi.phylo4<-reorder(ncbi.phylo4,order="preorder")
  #ncbi.phylo4.all.labels<-reorder(ncbi.phylo4.all.labels,order="preorder")
 # ncbi.phylo<-structure(list(edge=phylo.matrix.for.ape,tip.label=tips.name,root.edge=original.root,Nnode=dim(phylo.matrix.for.ape)[1]-length(tips.name)),class="phylo")
  #ncbi.phylo<-reorder(ncbi.phylo)
  if(download) {
    setwd("..")
    system("rm -r ncbi_taxdmp")
  }
  return(list(tree.phylo4=ncbi.phylo4,phylo.matrix=phylo.matrix))
}
