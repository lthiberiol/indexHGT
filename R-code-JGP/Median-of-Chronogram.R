################################################################################################
### Massachusetts Institute of Technology                                                      #
### Department of Earth, Atmospheric & Planetary Sciences, Fournier Lab                        #
### HGT-analysis-script                                                                        #
### Purpose: To produce a median date estimate from 9900 Trees from model BB for Erik Tamre    #
### Authors: Jack G. Payette, M.S. (https://orcid.org/0000-0001-6479-4557)                     #
### Authors2: L Thiberio Rangel, Ph.D. & Professor Gregory P. Fournier, Ph.D.                  #
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

### Source functions ###

#Function to load index files for use in HGT filter matching.
#REQUIRES SETTING WORKING DIRECTORY!
LoadHGTrun <- function(nameofrun2load){
  nameofrun<-nameofrun2load
  filename1 <- "PassedTreesIndex.txt"
  filename1 = paste(nameofrun,filename1)
  HGTRunIndex<-read.csv(file=filename1)
  assign(paste0("PassedTreesIndex",(as.character(nameofrun)),sep=""),HGTRunIndex,envir = .GlobalEnv)
};
# LoadHGTRun(1)

### LoadDateDistNode Function: Example underneath.
LoadDateDistNodeComments <-function(LabelDistFilePath,DateDistFilePath,UniqueModelCode,NodeNumberToFilter,Comments){
  ### Load Reference Label Tree as MyTreeLabels
  label_input = LabelDistFilePath
  MyTreeLabels = read.tree(label_input)
  class(MyTreeLabels) <- "multiPhylo"
  MyTreeLables <<- MyTreeLabels
  Nnodes <<- MyTreeLabels$Nnode
  internalNodeLabels <<- data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  Nnodes <- MyTreeLabels$Nnode
  print(paste0("MyTreeLabels object loaded from tree labels file with total Nnodes=",Nnodes))
  model <- read.tree(DateDistFilePath);
  class(model) <- "multiPhylo";
  MyTreeInput <- model
  assign(paste0("NumberOfTreesIn",UniqueModelCode),length(model),envir = .GlobalEnv)
  assign(paste0(UniqueModelCode),model,envir = .GlobalEnv);
  assign(paste0(UniqueModelCode,"DateDist",sep=""),DateDistFilePath,envir = .GlobalEnv);
  assign(paste0("NodeNumberToFilter"),NodeNumberToFilter,envir = .GlobalEnv);
  print("Vars: MyTreeLabels,Nnodes,internalNodeLabels,DateDist,ModelCode,NumTrees,Node#ToFilter")
  n=0
  modelsamplednodenum=0
  for (i in MyTreeInput)
  {
    n=n+1
    modelsamplednodenum[[n]] = as.numeric(i$node.label[with(internalNodeLabels,N[Node==(NodeNumberToFilter)])])
  };
  assign(paste0(UniqueModelCode,NodeNumberToFilter), modelsamplednodenum,envir = .GlobalEnv);
  print(paste0(paste0("Tree filtered on Node#",NodeNumberToFilter)," and vector created!"))
  print(paste0(paste0("Data Loaded for Model: ",UniqueModelCode,sep=" "),Comments));
}
# LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",177,"cyanobacteria crown")

### LoadDateDist Function: Example underneath.
LoadDateDist <-function(DateDistFilePath,UniqueModelCode){
  model <- read.tree(DateDistFilePath);
  class(model) <- "multiPhylo";
  assign(paste0(UniqueModelCode),model,envir = .GlobalEnv);
  assign(paste0(UniqueModelCode,"datedist",sep=""),DateDistFilePath,envir = .GlobalEnv);
}
# LoadDateDist("modeldata/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB")

### Analysis functions ###

### PrintNodeAge Function: New Version w/ MyTreeLabels 9-18-19 Example underneath.
PrintNodeAge <- function(MyTreeLabels,MyTreeInput,TreeNumber,StandardNodeNumber){
  Nnodes <- MyTreeLabels$Nnode
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeLabels$node.label)))
  NodeAge<-as.numeric(MyTreeInput[[TreeNumber]]$node.label[with(internalNodeLabels,N[Node==(StandardNodeNumber)])])
  returnValue(NodeAge);};
# PrintNodeAge(MyTreeLabels,BB,1,177) #Example use case

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

PlotHistModelNode<-function(UniqueModelCode,NodeNumberToFilter,ModelNode){
  #Can't get this to work yet: assign("ModelNode",as.name(paste0(UniqueModelCode,NodeNumberToFilter)),envir=.GlobalEnv)
  ModelTrees<-length(ModelNode)
  pdf(paste0(paste0(paste0(paste0("Histogram-model-",UniqueModelCode),"-Nodes-Age-"),NodeNumberToFilter),".pdf"), width=10, height=6)
  hist(ModelNode,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",
       ylab=paste0(paste0("Percent of Total (Frequency; N=",ModelTrees),")"),main=(paste0(paste0(paste0("Histogram of Age Estimates for Node#",NodeNumberToFilter)," for PhyloBayes model:\n"),UniqueModelCode)));
  dev.off();
}

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit");

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

####################################################################
#Select Trees that passed multiple HGT constraints:
# HGT[n] + HGT[n+1]   ### TO-DO recode this in a loop (maybe?)
# Get the runs you want to match...
#TO-DO: Recode / Figure out argument to pass in to match_df to avoid any warning ,on="x"??
#Existing script re-orders columns so what was the index column becomes the var.
#TO-DO: Recode this in Python using a better intersection algorithm!!!
### ANALYZE FOR INDIVIDUAL MODEL:

### BB long: 9,990 TREES
LoadDateDist("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/Cyano_modelBB_long_ugam_bd_7_20_sample.datedist","BB")
# "All Initial Trees.datedist"

#setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/")
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBB_long")
#setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/deliverables/RunBB_long-10-15-19")

################################################################################################
### DEFINE FUNCTIONS ###########################################################################
################################################################################################

# From original script to plot the median node age on a single tree from a given datedist file
# Writes to a table and PDF
PlotPassedTrees <- function(datedistfile)
{
  AllTrees <- read.tree(file=datedistfile)
  class(AllTrees) <- "multiPhylo"
  AllTreesDates <- AllTrees
  # Start
  conditionalDates = matrix(ncol = length(AllTreesDates[[1]]$node.label),nrow=9999)  #(optional) a vector of mode character giving the names of the nodes.
  conditionalEdges = matrix(ncol = length(AllTreesDates[[1]]$edge.length),nrow=9999)  #(optional) a numeric vector giving the lengths of the branches given by edge.
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
### DATA ANALYSIS            ###################################################################
################################################################################################

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/deliverables/Erik-10-24-19")

PlotPassedTrees("Cyano_modelBB_long_ugam_bd_7_20_sample.datedist");

#Outputs: PDF, table, and new tree file for Median Date Estimates from all trees/nodes in model.

##############################################
### END SCRIPT EXECUTION #####################
##############################################
