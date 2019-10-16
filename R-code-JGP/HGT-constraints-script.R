################################################################################################
### Massachusetts Institute of Technology                                                      #
### Department of Earth, Atmospheric & Planetary Sciences, Fournier Lab                        #
### HGT-constraints-script                                                                     #
### Purpose: Filtering/Post sampling a datedist file & tree(s) for given HGT contraints        #
### Authors: Jack G. Payette, M.S. (https://orcid.org/0000-0001-6479-4557)                     #
### Authors2: Cara Magnabosco, Ph.D. & Professor Gregory P. Fournier, Ph.D.                    #
################################################################################################

### Load Required Libraries
library(ape)
library(coda)
library(mcmcplots)
library(phangorn) #Tip for installing: restart R studio or use NO option on, from sources that need compilation!
library(plyr)
library(dplyr) #Diagnostic: Some function masking.
library(ggplot2)
library(reshape2)
library(tools) #Added for filtering
library(lattice) #Added for plotting
library(packrat) #This R studio project uses 'packrat' to manage packages independently/locally.

# Environment diagnostics
sessioninfo::session_info()
sessioninfo::package_info() #Diagnostic: Some of the packages might have a loaded and on-disk path mismatch.
#NOW Running on: R version 3.6.1 (2019-07-05) & ape 5.3

################################################################################################
### DEFINE FUNCTIONS ###########################################################################
################################################################################################
### Filter/Post sample a given datedist for ONLY Trees that pass HGT constraints
#Inputs: - multiPhylo datedist file = MyTree ***
#        - HGT constraint for internal node # ordered index
#        - a,d are passed in from HGTnode_input file which maps to tree labels file
# *** Ancestor node OLDER than Descendant node; a > d ***
#        a row in HGT index file = run name
#Outputs:a "PassedTrees.datedist" with prefix of run name.
#TO-DO COMPLETED FIX FOR ERROR ON <2,=1 TREES IN LOOP!!!
#TO-DO allow a project name to be passed into this function and written to output files
#TO-DO allow for different A>D constraints, beyond just internal node order, perhaps an option =1,2,3 for internal,standard,ranger nodes
################################################################################################

# FIX 10-8-19 LIMITED TO 9999 TREES CURRENTLY!
SampleTrees0 <- function(tree0,a,d,nameofrun,print=T){
  MyTree <- tree0
  n=0
  tmp = 0
  PassedTrees <- vector("list", 2) #create an empty list of 2 elements to store PassedTrees
  class(PassedTrees) <- "multiPhylo" #make this list a multiPhylo object
  nPassedTree <- vector("list", 1) #create an empty list of 2 elements to store PassedTrees
  class(nPassedTree) <- "phylo" #make this list a multiPhylo object
  write.tree(MyTree,"All Initial Trees.datedist")
  AllTrees <- read.tree(file="All Initial Trees.datedist")
  class(AllTrees) <- "multiPhylo"
  TreePassedConstraint<-as.vector(c(seq(AllTrees)))
  #Dates
  conditionalDates = matrix(ncol = length(AllTrees[[1]]$node.label),nrow=9999)  #(optional) a vector of mode character giving the names of the nodes.
  conditionalEdges = matrix(ncol = length(AllTrees[[1]]$edge.length),nrow=9999)  #(optional) a numeric vector giving the lengths of the branches given by edge.
  # For loop to evaluate one HGT
  nPassedTrees <-c(0)
  for (i in MyTree){
    tmp = tmp + 1
    dates = as.numeric(i$node.label)
    A = dates[a]
    D = dates[d]
    TreePassedConstraint[tmp]<-" "
    if (A>D){
      print(tmp)
      n=n+1
      nPassedTrees=nPassedTrees+1
      TreePassedConstraint[tmp]<-1 #Write to a text file with '1' for PassedTrees
      conditionalDates[n,]=dates #Turned OFF because no data manipulation just filtering
      conditionalEdges[n,]=i$edge.length
      PassedTrees[[n]]=i
      nPassedTree=i
    }
  }
  if (is.null(nPassedTree[[1]])==TRUE)
  {  nPassedTrees<-""
  filename1 <- " PassedTrees.datedist"
  filename1 = paste0(nameofrun,filename1,sep="")
  ZeroPassedTree = ";"
  write(ZeroPassedTree,file=filename1)
  filename2 <- " PassedTreesIndex.txt"
  filename2 = paste0(nameofrun,filename2,sep="")
  write.csv(TreePassedConstraint,file=filename2)
  return(nPassedTrees)

  }
  else if ( (nPassedTrees<2)==FALSE )
  {
    filename1 <- " PassedTrees.datedist"
    filename1 = paste0(nameofrun,filename1,sep="")
    write.tree(PassedTrees,file=filename1)
    filename2 <- " PassedTreesIndex.txt"
    filename2 = paste0(nameofrun,filename2,sep="")
    write.csv(TreePassedConstraint,file=filename2)
    return(PassedTrees)
  }
  else if ( (nPassedTrees=1)==TRUE )
  {
    filename1 <- " PassedTrees.datedist"
    filename1 = paste0(nameofrun,filename1,sep="")
    write.tree(nPassedTree,file=filename1)
    filename2 <- " PassedTreesIndex.txt"
    filename2 = paste0(nameofrun,filename2,sep="")
    write.csv(TreePassedConstraint,file=filename2)
    return(nPassedTree)
  }

};
SampleTreesStandard0 <- function(MyTree,Ancestor,Descendant,i){
  a<- with(internalNodeLabels,N[Node==(Ancestor)])
  d<- with(internalNodeLabels,N[Node==(Descendant)])
  SampleTrees0(MyTree,a,d,i) #Given a set of Trees and constraints, filter and write results to file
};
# Examples for BA
#SampleTreesStandard0(BA,292,268,1) #HGT1
#SampleTreesStandard0(BA,293,324,2) #HGT2

#SampleTreesStandard0(BA,292,268,1) #HGT1
#SampleTreesStandard0(BA,290,177,10) #HGT10
#SampleTreesStandard0(BA,290,179,11) #HGT11
#SampleTreesStandard0(BA,287,177,13) #HGT13
#SampleTreesStandard0(BA,287,179,14) #HGT14

#SampleTreesStandard0(BB,292,268,1) #HGT1
#SampleTreesStandard0(BB,290,177,10) #HGT10
#SampleTreesStandard0(BB,290,179,11) #HGT11
#SampleTreesStandard0(BB,287,177,13) #HGT13
#SampleTreesStandard0(BB,287,179,14) #HGT14

### LoadDateDist Function: Example underneath.
LoadDateDist <-function(DateDistFilePath,UniqueModelCode){
  model <- read.tree(DateDistFilePath);
  class(model) <- "multiPhylo";
  assign(paste0(UniqueModelCode),model,envir = .GlobalEnv);
  assign(paste0(UniqueModelCode,"datedist",sep=""),DateDistFilePath,envir = .GlobalEnv);
};
#LoadDateDist("modeldata/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB")

### Loop over all HGT constraints and for each, save PassedTrees in a new datedist file with prefix# as row of HGT constraint.
#FIXED ON 10-8-19: Have better failure/error handling for returning 0 or 1 trees!!!
#NOTE TO-DO: Include the Node naming checks within this function!
#LOOP FOR each i in HGTnode
#Inputs: a multiPhylo datedist file = MyTree = tree0
#        a HGT index
#        a ProjectName = nameofrun
#Outputs: to console, table/Text
#Outputs: DATA from sub function called SampleTrees() outputs to datedist file PassedTrees.txt
FilterTreeHGT <- function(tree2,hgts,nameofrun,print=T){
  #Initialize w/ variables
  MyTree <- tree2
  HGT<-hgts
  #PassedHGT <- data.frame()
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  n=0
  tmp = 0
  PassedHGT <- vector("double")
  for (i in HGT)
  {
    tmp = tmp + 1
    n=n+1
    #print(tmp) #turn off since print=T and sub function prints the results, i.e. tree looped over
    #given Ancestor value find Descendant value in HGT table
    Ancestor<-HGTnode[i,1]
    Descendant<-HGTnode[i,2]
    a<- with(internalNodeLabels,N[Node==(Ancestor)])
    d<- with(internalNodeLabels,N[Node==(Descendant)])
    SampleTrees0(MyTree,a,d,i) #Given a set of Trees and constraints, filter and write results to file
    print(paste(i,c("HGT constraint run complete!"),sep=" "))
    PassedHGT[[n]]=  ReadPassedTrees(i) #Error: Error in if (all(phy$node.label == "")) phy$node.label <- NULL :
    #missing value where TRUE/FALSE needed Not the cleanest results function/table but can adjust.
  }
  write.table(PassedHGT,paste(nameofrun,"PassedHGTResults.txt"))
  return(PassedHGT) #Need to fix since this returns blank!
};

# Read a given datedist file from an HGT filter/sample run and compare to original datedist and output
# Input: a run number corresponding to HGT index (PassedTree.txt index files must exist on disk and do have spaces)
# Outputs: text file
ReadPassedTrees <- function(nameofrun,print=T){
  AllTrees <- read.tree(file="All Initial Trees.datedist")
  class(AllTrees) <- "multiPhylo"
  #Create new function to load back trees from output to test
  #filename1 <- " PassedTrees.datedist"
  #filename1 = paste0(nameofrun,filename1,sep="")
  #PassedTrees1 <- read.tree(file=filename1)
  filename1 <- " PassedTreesIndex.txt"
  filename1 = paste0(nameofrun,filename1,sep="")
  PassedTrees <- read.csv(filename1)
  #PassedTrees <- read.csv("1 PassedTreesIndex.txt")
  #PassedTrees <- read.csv("10 PassedTreesIndex.txt")
  #PassedTrees <- read.csv("2 PassedTreesIndex.txt")
  TreeThatPassed <- data.frame(PassedTrees$x)
  NumberOfPassedTrees <- as.numeric(lengths(filter(TreeThatPassed, PassedTrees.x == "1")))
  if (is.na(as.numeric(NumberOfPassedTrees))==FALSE){Number <- NumberOfPassedTrees}
  else {Number <-c(0)}
  checkoutput <-paste(paste("Results:",paste(Number,length(AllTrees),sep=" of ")),"Trees passed HGT constraints specified from datedist & file inputs.")
  return(checkoutput)
};
#Examples for BA
#ReadPassedTrees(1) #expect 1
#ReadPassedTrees(10) #expect zero
#ReadPassedTrees(2) #expect 10

### Check & Fail for mis match in number of internal nodes between trees in original datedist and from Reference Label Tree
#NOTE TO-DO: Have a way of checking order not just node number
#Inputs: a multiPhylo datedist file = MyTree = tree1 | a Reference Label Tree from file = MyTreeLabels = tree2
#Output: Returns NULL if Passed, Returns Failed Tree i from datedist if Assumption not met
CheckFailInternalNodeN <- function(tree1,tree2,print=T){
  tmp =0
  FailedTreei=vector("list")
  for (i in tree1) {
    tree0 <- i
    tmp = tmp+1
    if (tree0$Nnode != tree2$Nnode)
    {FailedNnode=tree0$Nnode
    FailedTreei<-c(rbind(FailedTreei,tmp))
    checkoutput<-cbind(FailedTreei,FailedNnode)
    return(checkoutput)
    print(i)
    }
  }
}; #Given tree 2 is reference label

# Check a given datedist file from an HGT filter/sample run
# Inputs: a run number corresponding to HGT index | a Refernce Label Tree = MyTreeLabels = tree2
CheckPassedTrees <- function(nameofrun,tree2,print=T){
  filename1 <- "PassedTrees.datedist"
  filename1 = paste(nameofrun,filename1)
  PassedTrees <- read.tree(file=filename1)
  class(PassedTrees) <- "multiPhylo"
  Cresult <- CheckFailInternalNodeN(PassedTrees,tree2)
  checkoutput <- is.null(Cresult)
  return(checkoutput)
};

#Function to load index files for use in HGT filter matching.
LoadHGTrun <- function(nameofrun2load){
  nameofrun<-nameofrun2load
  filename1 <- "PassedTreesIndex.txt"
  filename1 = paste(nameofrun,filename1)
  HGTRunIndex<-read.csv(file=filename1)
  assign(paste0("PassedTreesIndex",(as.character(nameofrun)),sep=""),HGTRunIndex,envir = .GlobalEnv)
};
#LoadHGTRun(1)

# HISTOGRAM PLOTTING AND FUCTION DEFINITION
histnode <- function(node2plot,filename0){
  PassedTreesHistogram<- read.tree(file=filename0)
  class(PassedTreesHistogram) <- "multiPhylo"
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  n=0
  nodelab=0
  for (i in PassedTreesHistogram)
  {
    n=n+1
    nodelab[[n]] = as.numeric(i$node.label[with(internalNodeLabels,N[Node==(node2plot)])])
  };
  filename3 <- paste(filename0,node2plot)
  write.csv(nodelab,file=paste("histogramofnode",filename3))
  pdf(paste(filename0,"HistogramTreeDates.pdf"), 7, 5)
  hist.default(nodelab,xlab="Node Label as Date in Ma",
               ylab=filename0,
               xlim=rev(c(2200,3600)))
  dev.off()
  return(hist.default(nodelab,xlab="Node Label as Date in Ma",
                      ylab=filename0,
                      xlim=rev(c(2200,3600))))

};
#histnode(177,"All Initial Trees.datedist")
#histnode(177,"1 PassedTrees.datedist")
#histnode(177,"4 PassedTrees.datedist") #Can't be executed for runs 1,3,5 which return NULL for HGT#4
#histnode(177,"All PassedTrees.datedist")

histnodedensity <- function(node2plot,filename0){
  PassedTreesHistogram<- read.tree(file=filename0)
  class(PassedTreesHistogram) <- "multiPhylo"
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  n=0
  nodelab=0
  for (i in PassedTreesHistogram)
  {
    n=n+1
    nodelab[[n]] = as.numeric(i$node.label[with(internalNodeLabels,N[Node==(node2plot)])])
  };
  filename3 <- paste(filename0,node2plot)
  write.csv(nodelab,file=paste("histogramofnode",filename3))
  pdf(paste(filename0,"DensityTreeDates.pdf"), 7, 5)
  plot(density(nodelab),main="Density of Node Age in Ma",ylab=filename0,xlim=rev(c(2200,3600)))
  dev.off()
  return(plot(density(nodelab),main="Density of Node Age in Ma",ylab=filename0,xlim=rev(c(2200,3600))))
};
#histnodedensity(177,"All Initial Trees.datedist")
#histnodedensity(177,"1 PassedTrees.datedist")
#histnodedensity(177,"4 PassedTrees.datedist") #Can't be executed for runs 1,3,5 which return NULL for HGT#4
#histnodedensity(177,"All PassedTrees.datedist")

################################################################################################
### IN DEV FUNCTIONS ###########################################################################
#Utility Fn # Requires separate MyTreeLabels call

### Function to Re-label a multiPhylo: MyTreeInput$tip.label <<- MyTreeLabelsNewTips
ReMapTipLabels <- function(MyTreeInput,MyTreeLabelsNewTips,UniqueModelCode){
  n=0
  class(MyTreeInput) <- "multiPhylo"
  MyTreeOutput=vector(mode="list",length=c(5));
  class(MyTreeOutput) <- "multiPhylo";
  for (i in MyTreeInput[]){
    n=n+1
    class(i) <- "phylo"
    class(MyTreeLabelsNewTips) <-"phylo"
    MyNewTipLabels <<- MyTreeLabelsNewTips$tip.label
    i$tip.label=MyTreeLabelsNewTips$tip.label
    MyTreeOutput[[n]]= i
    class(MyTreeOutput) <- "multiPhylo"
  }
  assign(paste0(UniqueModelCode),MyTreeOutput,envir = .GlobalEnv)
};
# ReMapTipLabels(BB,MyTreeLabelsNewTips,"BBnewtips")

### PrintNodeAge Function: New Version w/ MyTreeLabels 9-18-19 Example underneath.
PrintNodeAge <- function(MyTreeLabels,MyTreeInput,TreeNumber,StandardNodeNumber){
  Nnodes <- MyTreeLabels$Nnode
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  NodeAge<-as.numeric(MyTreeInput[[TreeNumber]]$node.label[with(internalNodeLabels,N[Node==(StandardNodeNumber)])])
  returnValue(NodeAge);};
# PrintNodeAge(MyTreeLabels,BB,1,177) #Example use case

### PrintRANGERNodeAge: New Version w/ Ranger Labels Remapping! Example underneath.
PrintRangerNodeAge <- function(MyTreeLabels,MyTreeLabelsRanger,MyTreeInput,TreeNumber,RangerNode){
  NnodesRanger <- MyTreeLabelsRanger$Nnode
  internalNodeLabelsRanger = data.frame(N = c(seq(1:NnodesRanger)),NodeNameRanger=c((MyTreeLabelsRanger$node.label)))
  Nnodes <- MyTreeLabels$Nnode
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  NodeAge<-as.numeric(MyTreeInput[[TreeNumber]]$node.label[with(internalNodeLabelsRanger,N[NodeNameRanger==(RangerNode)])])
  returnValue(NodeAge)};
# PrintRangerNodeAge(MyTreeLabels,MyTreeLabelsRanger,BB,1,'n61')

# From original script to plot the median node age on a single tree from a given datedist file
# Writes to a table and PDF
PlotPassedTrees <- function(datedistfile){
  AllTrees <- read.tree(file=datedistfile)
  class(AllTrees) <- "multiPhylo"
  AllTreesDates <- AllTrees
  # Start
  conditionalDates = matrix(ncol = length(AllTreesDates[[1]]$node.label),nrow=3000)  #(optional) a vector of mode character giving the names of the nodes.
  conditionalEdges = matrix(ncol = length(AllTreesDates[[1]]$edge.length),nrow=3000)  #(optional) a numeric vector giving the lengths of the branches given by edge.
  n=0
  tmp = 0
  foo <- vector("list", 2) #create an empty list of 2 elements
  class(foo) <- "multiPhylo" #make this list a multiPhylo object
  treenums <-vector("list")
  for (i in AllTreesDates){
    tmp = tmp + 1
    dates = as.numeric(i$node.label)
    print(tmp)
    n=n+1
    conditionalDates[n,]=dates
    conditionalEdges[n,]=i$edge.length
    foo[[n]]=i
  }
  ## Matrix containing node ages (columns) of trees (rows) that meet the conditional constraint.
  conditionalDates=conditionalDates[complete.cases(conditionalDates),]
  conditionalEdges=conditionalEdges[complete.cases(conditionalEdges),]
  # Calculations for Dates
  medPost = apply(conditionalDates,2,FUN = median)
  medEdge = apply(conditionalEdges,2,FUN=median)
  meanEdge = colMeans(conditionalEdges)
  maxPost = apply(conditionalDates,2,FUN = max)
  minPost = apply(conditionalDates,2,FUN = min)
  ciPost = paste(maxPost,minPost, sep=" - ")
  # Write dates to file
  outTable = t(rbind(medPost,ciPost))
  write.csv(outTable,file="All Passed MedianTrees Output table.txt")
  # Add dates to new tree & node label for Median Posterior Age
  newTree = AllTreesDates[[1]] # QUESTION: does this need to be manually set to the correct trees outputted by the previous loop? I think so...
  newTree$node.label=medPost
  #IF CONFIDENCE INTERVAL DESIRED USE THIS OPTION COMMENTED OUT
  #newTree$node.label=paste(apply(conditionalDates,2,FUN=max),apply(conditionalDates,2,FUN=min),sep=" - ")
  # Adjust edges ? TO-DO: Understand this better???
  tempTree = AllTreesDates[[1]] #Same question as above
  tempTree$node.label = seq(length(tempTree$tip.label)+1,length(tempTree$tip.label)+length(tempTree$tip.label)-1 )
  v = rep(0,length(tempTree$edge.length))
  # Adjust edge and node labels ???
  edge = tempTree$edge
  nodelabels = tempTree$node.label
  for (r in 1:nrow(edge)){
    outerNodeIndex = edge[r,1]   #outerNodeIndex = edge[1,1]
    outerNode = newTree$node.label[which(nodelabels==outerNodeIndex)]
    if (edge[r,2]>length(tempTree$tip.label)){
      innerNodeIndex = edge[r,2]
      innerNode = newTree$node.label[which(nodelabels==innerNodeIndex)]
    }
    else{
      innerNode = 0
    }
    v[r] = outerNode-innerNode
  }
  #Store edge length in a variable
  newTree$edge.length = v

  #To-Do: Rewrite keys from file
  #DEBUG OFF: taxonomyKey = read.delim("phototrophy_key_no_bin_key.txt",header=FALSE)
  tip = newTree$tip.label
  #DEBUG OFF: newTree$tip.label=mapvalues(tip, from=taxonomyKey$V3, to=as.character(taxonomyKey$V4))

  #TO-DO expand upon this!

  #Begin PDF
  pdf("AllMedianTreeDates.pdf", 7, 5)

  #TO-DO: Adjust these plot methods and understand them better!

  #plot1<- plot(newTree,show.node.label = FALSE,cex = 0.5,show.tip.label = TRUE)
  plot(newTree,show.node.label = TRUE,cex = 0.5,show.tip.label = TRUE)
  #axisPhylo()
  #plot1<- plot(newTree,show.node.label = FALSE,cex = 0.5,show.tip.label = TRUE)

  #tmp = newTree
  #tmp$node.label=seq(tmp$node.label)
  #plot(tmp,show.node.label = TRUE,cex = 0.5,show.tip.label = FALSE)
  #plot2<- plot(tmp,show.node.label = TRUE,cex = 0.5,show.tip.label = FALSE)
  axisPhylo()

  dev.off()
  write.tree(newTree, file="AllPassedMedianTreeDates.tree")
  return("AllMedianTreeDates.pdf")
};

################################################################################################
### BEGIN SCRIPT EXECUTION #####################################################################
################################################################################################
### Set Working Directory from MIT Dropbox ###
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit");
#setwd("/Users/Jack Payette/Dropbox (MIT)/R-code-mit/R-code-mit/")
print(getwd());

################################################################################################
###### LOAD NODE MAPPING AND HGT INDEX DATA FILES
###### TO-DO: Import Ranger node names Tree data to use as Reference Label Tree (with other HGT master file)
### Load Reference Label Tree as MyTreeLabels #Fixed class on 9-27-19 to be phylo NOT multiPhylo
label_input = "modeldata/Cyano_modelBB_ugam_bd_7_20_sample.labels"
MyTreeLabels = read.tree(label_input); class(MyTreeLabels) <- "phylo";
Nnodes <- MyTreeLabels$Nnode
internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
print(Nnodes); print("MyTreeLabels object loaded from tree labels file!");
### Load Ranger Label Tree as MyTreeLabelsRanger #Fixed class on 9-27-19 to be phylo NOT multiPhylo
label_input_ranger = "modeldata/RANGER-species-tree-named-nodes"
MyTreeLabelsRanger = read.tree(label_input_ranger);class(MyTreeLabelsRanger) <- "phylo";
NnodesRanger <- MyTreeLabelsRanger$Nnode
internalNodeLabelsRanger = data.frame(N = c(seq(1:NnodesRanger)),NodeNameRanger=c((MyTreeLabelsRanger$node.label)))
print(NnodesRanger); print("MyTreeLabelsRanger object loaded from tree labels file!");
### Re-label / Map Tips as MyTreeLabelsNewTips from MyTreeLabels
MyTreeLabelsNewTips <- MyTreeLabels; class(MyTreeLabelsNewTips) <- "phylo";
taxonomyKey = read.csv("modeldata/Kelsey_CYPL_key_payette.csv",header=TRUE,sep=",")
tip = MyTreeLabels$tip.label
MyTreeLabelsNewTips$tip.label=mapvalues(tip, from=taxonomyKey$ShortName, to=as.character(taxonomyKey$Code))
## Set these to MyTreeLabels & MyTreeLablesRanger
# MyTreeLabels$tip.label <- MyTreeLabelsNewTips$tip.label
# MyTreeLabelsRanger$tip.label <- MyTreeLabelsNewTips$tip.label

### HGT index as data frame and HGT atomic vector index
####### TO-DO: Debug error in read.table, which will occur because of file format.
HGTnode_input <- read.csv(file="modeldata/HGTindex_highconfidence_10_2_19.txt", header=TRUE, sep=",",fill=TRUE)
# N.B. Carriage return/mac vs PC text file can throw error, open file, reformat, try again
# HGTnode index and atomic vector HGT also created
HGTnode <- data.matrix(HGTnode_input[,])
HGT <- as.vector(seq(1,length(HGTnode_input[,1])),mode="numeric")
print("HGT node index loaded!")
print(HGTnode_input)
############# Build Node Labels list from Reference Label Tree: MyTreeLabels
#TO-DO: Re-code with index from RANGER internal node labelling ##### TO-DO save this file as an index
# Number of Nodes created from Reference Label Tree: MyTreeLabels
Nnodes <- MyTreeLabels$Nnode
internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
print(Nnodes)

################################################################################################
###### LOAD DATA FILES from datedist 'as' multiPhylo (ape) ######
#TO-DO COMPLETED: Better manage file inputs from a function - have project/run names passed through to output files
#datedist_input = "/data/" #run
#LoadDateDist("r1_1.2akinete_sample.datedist","run0")
#LoadDateDist("data/Cyano_modelBC_ugam_bd_7_20_sample.datedist","run1")
#LoadDateDist("data/Cyano_modelBE_ugam_bd_7_20_sample.datedist","run2")
#LoadDateDist("data/Cyano_modelBE_ugam_nobd_7_20_sample.datedist","run3")
#LoadDateDist("data/Cyano_modelBE_ugam_bd_longrun_7_20_sample.datedist","run4")
#LoadDateDist("data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","run5")
#LoadDateDist("data/Cyano_modelBG_ugam_bd_7_20_sample.datedist.datedist","run6")
#
LoadDateDist("modeldata/Cyano_modelBA_ugam_bd_7_20_sample.datedist","BA")
LoadDateDist("modeldata/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB")
LoadDateDist("modeldata/Cyano_modelBC_ugam_bd_7_20_sample.datedist","BC")
LoadDateDist("modeldata/Cyano_modelBD_ugam_bd_7_20_sample.datedist","BD")
LoadDateDist("modeldata/Cyano_modelBE_ugam_bd_7_20_sample.datedist","BE")
LoadDateDist("modeldata/Cyano_modelBF_ugam_bd_7_20_sample.datedist","BF")
LoadDateDist("modeldata/Cyano_modelBG_ugam_bd_7_20_sample.datedist","BG")
#LoadDateDist(" "," ")

################################################################################################
###### CHECKSUM ###################################################
# Check Assumption that Internal Node # == for all Trees compared to MyTreeLabels
CheckFailInternalNodeN(MyTree,MyTreeLabels) # Expect to RETURN NULL
MyTree[[1]]$Nnode == MyTreeLabels$Nnode # Expect to RETURN TRUE Validate Internal node numbers to ensure match
print(Nnodes) #Print N expected nodes

CheckFailInternalNodeN(BA,MyTreeLabels)
CheckFailInternalNodeN(BB,MyTreeLabels)
CheckFailInternalNodeN(BC,MyTreeLabels)
CheckFailInternalNodeN(BD,MyTreeLabels)
CheckFailInternalNodeN(BE,MyTreeLabels)
CheckFailInternalNodeN(BF,MyTreeLabels)
CheckFailInternalNodeN(BG,MyTreeLabels)

BA[[1]]$Nnode == MyTreeLabels$Nnode
BB[[1]]$Nnode == MyTreeLabels$Nnode
BC[[1]]$Nnode == MyTreeLabels$Nnode
BD[[1]]$Nnode == MyTreeLabels$Nnode
BE[[1]]$Nnode == MyTreeLabels$Nnode
BF[[1]]$Nnode == MyTreeLabels$Nnode
BG[[1]]$Nnode == MyTreeLabels$Nnode

################################################################################################

### EXECUTE LOOP OVER ALL HGT AND FILTER FOR TREES THAT MEET CONSTRAINT

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBA/")
FilterTreeHGT(BA,HGT,"BA_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBB/")
FilterTreeHGT(BB,HGT,"BB_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBC/")
FilterTreeHGT(BC,HGT,"BC_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBD/")
FilterTreeHGT(BD,HGT,"BD_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBE/")
FilterTreeHGT(BE,HGT,"BE_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBF/")
FilterTreeHGT(BF,HGT,"BF_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBG/")
FilterTreeHGT(BG,HGT,"BG_model_cyanoclock_highpass")

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/RunBE_long/")
LoadDateDist("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/data/Cyano_modelBE_ugam_bd_longrun_7_20_sample.datedist","BElong")
FilterTreeHGT(BElong,HGT,"BE_longrun_model_cyanoclock_highpass") #5700 trees

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBE_nobd/")
LoadDateDist("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/data/Cyano_modelBE_ugam_nobd_7_20_sample.datedist","BEnobd")
FilterTreeHGT(BEnobd,HGT,"BE_nobd_model_cyanoclock_highpass") #1400 trees

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBE_ugambd")
LoadDateDist("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/data/Cyano_modelBE_ugam_bd_7_20_sample.datedist","BEugambd")
FilterTreeHGT(BEugambd,HGT,"BE_ugambd_model_cyanoclock_highpass") #1500 trees

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBB_long")
LoadDateDist("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/Cyano_modelBB_long_ugam_bd_7_20_sample.datedist","BBlong")
FilterTreeHGT(BBlong,HGT,"BB_long_ugam_bd_model_cyanoclock_highpass") #9900 trees

### CHECK A GIVEN HGT CONSTRAINT THAT IT HAS SAME NUMBER INTERNAL NODES
# TO-DO RECODE AS LOOP FOR EACH I IN HGT:

CheckPassedTrees(1,MyTreeLabels)
CheckPassedTrees(2,MyTreeLabels)
CheckPassedTrees(3,MyTreeLabels)
CheckPassedTrees(4,MyTreeLabels)
CheckPassedTrees(5,MyTreeLabels)
CheckPassedTrees(6,MyTreeLabels)
CheckPassedTrees(7,MyTreeLabels)
CheckPassedTrees(8,MyTreeLabels)
CheckPassedTrees(9,MyTreeLabels)
CheckPassedTrees(10,MyTreeLabels)
CheckPassedTrees(11,MyTreeLabels)
CheckPassedTrees(12,MyTreeLabels)
CheckPassedTrees(13,MyTreeLabels)
CheckPassedTrees(14,MyTreeLabels)
CheckPassedTrees(15,MyTreeLabels)
CheckPassedTrees(16,MyTreeLabels)

histnode(177,"1 PassedTrees.datedist")
histnode(177,"2 PassedTrees.datedist")
histnode(177,"3 PassedTrees.datedist")
histnode(177,"4 PassedTrees.datedist")
histnode(177,"5 PassedTrees.datedist")
histnode(177,"6 PassedTrees.datedist")
histnode(177,"7 PassedTrees.datedist")
histnode(177,"8 PassedTrees.datedist")
histnode(177,"9 PassedTrees.datedist")
histnode(177,"10 PassedTrees.datedist")
histnode(177,"11 PassedTrees.datedist")
histnode(177,"12 PassedTrees.datedist")
histnode(177,"13 PassedTrees.datedist")
histnode(177,"14 PassedTrees.datedist")
histnode(177,"15 PassedTrees.datedist")
histnode(177,"16 PassedTrees.datedist")

histnodedensity(177,"1 PassedTrees.datedist")
histnodedensity(177,"2 PassedTrees.datedist")
histnodedensity(177,"3 PassedTrees.datedist")
histnodedensity(177,"4 PassedTrees.datedist")
histnodedensity(177,"5 PassedTrees.datedist")
histnodedensity(177,"6 PassedTrees.datedist")
histnodedensity(177,"7 PassedTrees.datedist")
histnodedensity(177,"8 PassedTrees.datedist")
histnodedensity(177,"9 PassedTrees.datedist")
histnodedensity(177,"10 PassedTrees.datedist")
histnodedensity(177,"11 PassedTrees.datedist")
histnodedensity(177,"12 PassedTrees.datedist")
histnodedensity(177,"13 PassedTrees.datedist")
histnodedensity(177,"14 PassedTrees.datedist")
histnodedensity(177,"15 PassedTrees.datedist")
histnodedensity(177,"16 PassedTrees.datedist")

##############################################
### END SCRIPT EXECUTION #####################
##############################################

####################################################################
### DRAFT CODE FOR ANALYSIS FOR INDIVIDUAL MODEL:
### BB_long / BEugambd
####################################################################

#Select Trees that passed multiple HGT constraints:
# HGT[n] + HGT[n+1]   ### TO-DO recode this in a loop (maybe?)
#Get the runs you want to match...
#TO-DO: Recode / Figure out argument to pass in to match_df to avoid any warning ,on="x"??
#Existing script re-orders columns so what was the index column becomes the var.
#TO-DO: Recode this in Python using a better intersection algorithm!!!

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/")
# NOTE! This must contain index files for the model/datedist you're working with!

for (i in HGT){LoadHGTrun(i)}

#HGT Group 0 All Initial Trees
#HGT Group 1: 4,8,15,7,16,9,3,6,12,5
#HGT Group 2:                       ,2
#HGT Group 3:                        ,14
#HGT Group 4:                           ,13
#HGT Group 5:                              ,1
#HGT Group 6:                                ,11
#HGT Group 7:                                  ,10

#1 {4,8,15,7,16,9,3,6,12,5}

### Match HGT 3,5,7 to see which trees pass all these constraints!
MatchAllPassed <- data.frame(match_df(match_df(PassedTreesIndex3,PassedTreesIndex5),PassedTreesIndex7,on=NULL)) #use reduce/map to make this work!
#TO-DO Fix this sloppy code and use reduce/map to make this work!

# Intersect HGT Group 1
MatchAllPassed1 <- data.frame(
  match_df(
    match_df(
      match_df(
        match_df(
          match_df(
            match_df(
              match_df(
                match_df(
                  match_df(PassedTreesIndex4,PassedTreesIndex8,on=NULL
                  ),PassedTreesIndex15,on=NULL
                ),PassedTreesIndex7,on=NULL
              ),PassedTreesIndex16,on=NULL
            ),PassedTreesIndex9,on=NULL
          ),PassedTreesIndex3,on=NULL
        ),PassedTreesIndex6,on=NULL
      ),PassedTreesIndex12,on=NULL
    ),PassedTreesIndex5,on=NULL)
);
# Check how many trees pass
length(MatchAllPassed1$x)
length((filter(MatchAllPassed1, x == "1")$x))
# Data manipulation
AllPassed <- MatchAllPassed1$X[ which(MatchAllPassed1$x %in% 1) ]
LoadDateDist("All Initial Trees.datedist","MyInitialTrees")
AllPassedTrees <- MyInitialTrees[c(AllPassed)]
# Print Passed Trees
print(AllPassedTrees)
print(AllPassed)
# Write Passed Trees down to disk
write.tree(AllPassedTrees,"Match HGT 1 PassedTrees.datedist")
write.csv(AllPassed,file="Match HGT 1 PassedTrees.txt")

# Load back results
LoadDateDist("Match HGT 1 PassedTrees.datedist","MatchAllTrees1")

# Print results
PrintNodeAge(MyTreeLabels,MatchAllTrees1,1,177)

histnode(177,"Match HGT 1 PassedTrees.datedist")
histnodedensity(177,"Match HGT 1 PassedTrees.datedist")

histnode(177,"All Initial Trees.datedist")
histnodedensity(177,"All Initial Trees.datedist")

##############################################
### SAVED RESULTS FOR DEBUGGING ##############
##############################################

### END SCRIPT RESULTS & NOTES
##############################################
