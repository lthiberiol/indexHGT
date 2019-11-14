################################################################################################
### Massachusetts Institute of Technology                                                      #
### Department of Earth, Atmospheric & Planetary Sciences, Fournier Lab                        #
### HGT-analysis-script                                                                        #
### Purpose: Analyzing filtered/Post sampled model datedist data for given HGT contraints      #
### Authors: Jack G. Payette, M.S. (https://orcid.org/0000-0001-6479-4557)                     #
### Authors2: L Thiberio Rangel, Ph.D. & Professor Gregory P. Fournier, Ph.D.                   #
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
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager");
BiocManager::install("qvalue");
install.packages('jaccard');
library('jaccard');

library(heatmap.plus);
#
#install.packages("IntClust")
#library(IntClust)

# Environment diagnostics
sessioninfo::session_info()
sessioninfo::package_info() #Diagnostic: Some of the packages might have a loaded and on-disk path mismatch.

### Source functions ###

#Function to load index files for use in HGT filter matching.
LoadHGTrun <- function(nameofrun2load){
  nameofrun<-nameofrun2load
  filename1 <- "PassedTreesIndex.txt"
  filename1 = paste(nameofrun,filename1)
  HGTRunIndex<-read.csv(file=filename1)
  assign(paste0("PassedTreesIndex",(as.character(nameofrun)),sep=""),HGTRunIndex,envir = .GlobalEnv)
};
  # LoadHGTRun(1)

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

# Function to build vectors from HGT index files to use for clustering analysis
LoadHGTrunVector<- function(nameofrun2load){
  nameofrun<-nameofrun2load
  filename1 <- "PassedTreesIndex.txt"
  filename1 = paste(nameofrun,filename1)
  HGTRunIndex<-read.csv(file=filename1)
  hgtrun <- as.vector(HGTRunIndex$x)
  hgtrun[is.na(hgtrun)] <- 0
  assign(paste0("hgt_",(as.character(nameofrun)),sep=""),hgtrun,envir = .GlobalEnv)
  assign(paste0("PassedTreesIndex",(as.character(nameofrun)),sep=""),HGTRunIndex,envir = .GlobalEnv)
}
  #for (i in HGT){LoadHGTrunVector(i)};

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

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/")
#setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBE_long")

###### BEGIN BB LONG ANALYSIS

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunBB_long")
#setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/deliverables/RunBB_long-10-15-19")

#Load data from Working Directory:
LoadDateDist("All Initial Trees.datedist","MyInitialTrees")
#Load HGT run index -- which trees passed or not from the above original file
for (i in HGT){LoadHGTrun(i)};

#Manual code for Jaccard Index
hgt1<-as.vector(PassedTreesIndex1$x)
hgt1[is.na(hgt1)] <- 0
hgt2<-as.vector(PassedTreesIndex2$x)
hgt2[is.na(hgt2)] <- 0
#Manual code for Jaccard Distance
hgt_1_2_j_dist <- ((1-jaccard(hgt1,hgt2,center=FALSE,px=NULL,py=NULL))*100)

#Code for melting before converting to vector!
#HGTpassMelt <- melt(HGTpass, varnames=c(hgt_1,hgt_2,hgt_3,hgt_4,hgt_5,hgt_6,hgt_7,hgt_8,hgt_9,hgt_10,hgt_11,hgt_12,hgt_13,hgt_14,hgt_15,hgt_16))

# Ingest vectors for each HGT index and build a new data table/matrix with trees (rows) and a binary 0|1 indicating pass or not for a given HGT constraint (columns)
for (i in HGT){LoadHGTrunVector(i)};
# Build data frame
HGTpass <- cbind(hgt_1,hgt_2,hgt_3,hgt_4,hgt_5,hgt_6,hgt_7,hgt_8,hgt_9,hgt_10,hgt_11,hgt_12,hgt_13,hgt_14,hgt_15,hgt_16);
HGTpass <- data.frame(HGTpass);
colnames(HGTpass)[1:16] <- c(1:16) #numeric column names for each HGT constraint
#write this new table down to disk
write.csv2(HGTpass,"RunBB_long_HGTpass.csv");
#build matrix
HGTpassM <-as.matrix(HGTpass)

#Produce heatmap!
heatmap <- heatmap.plus(HGTpassM)

#example: lapply(seq_along(xs), function(i) {})
#doesn't work: lapply(c(colnames(HGTpass),colnames(HGTpass)), jaccard.test.exact)
#doesn't work: lapply(HGTpass[,x],function(i){jaccard(HGTpass[,i],HGTpass[,(i+1)])})
#example: lapply(names(HGTpass),function(nm){jaccard(HGTpass[nm],HGTpass[(nm)])});

hgt1ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt2ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[2],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt3ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[3],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt4ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[4],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt5ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[5],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt6ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[6],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt7ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[7],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt8ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[8],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt9ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[9],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt10ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[10],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt11ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[11],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt12ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[12],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt13ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[13],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt14ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[14],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt15ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[15],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
hgt16ji <- lapply(seq_along(HGTpass),function(i){ jaccard(HGTpass[16],HGTpass[i],center=FALSE,px=NULL,py=NULL) } );
HGTpassJi <- cbind(hgt1ji,hgt2ji,hgt3ji,hgt4ji,hgt5ji,hgt6ji,hgt7ji,hgt8ji,hgt9ji,hgt10ji,hgt11ji,hgt12ji,hgt13ji,hgt14ji,hgt15ji,hgt16ji);
write.csv2(HGTpassJi,file="JaccardIndexHGTpassBBlong.csv")


#changing Jaccard Index: hgt1j <- lapply(seq_along(HGTpass),function(i){jaccard(HGTpass[1],HGTpass[i])});
#to Jaccard Distance (dis-similar): hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt2j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[2],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt3j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[3],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt4j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[4],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt5j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[5],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt6j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[6],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt7j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[7],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt8j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[8],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt9j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[9],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt10j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[10],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt11j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[11],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt12j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[12],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt13j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[13],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt14j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[14],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt15j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[15],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt16j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[16],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );

#HGTpassJ <- cbind(hgt1j,hgt2j,hgt3j,hgt4j,hgt5j,hgt6j,hgt7j,hgt8j,hgt9j,hgt10j,hgt11j,hgt12j,hgt13j,hgt14j,hgt15j,hgt16j);
#write.csv2(HGTpassJ,file="JaccardIndexHGTpassBBlong.csv")

HGTpassJd <- cbind(hgt1j,hgt2j,hgt3j,hgt4j,hgt5j,hgt6j,hgt7j,hgt8j,hgt9j,hgt10j,hgt11j,hgt12j,hgt13j,hgt14j,hgt15j,hgt16j);
write.csv2(HGTpassJd,file="JaccardDistanceHGTpassBBlong.csv")
# Other code below to build a different data object of the same Jaccard Distance data
HGTpassJaccard <- data.frame(HGTpassJd)
colnames(HGTpassJaccard)[1:16] <- c(1:16)
HGTpassJaccard <- matrix(HGTpassJaccard)

#Key code to create matrix and heatmap
HGTpassJaccardDistance <- cbind(hgt1j,hgt2j,hgt3j,hgt4j,hgt5j,hgt6j,hgt7j,hgt8j,hgt9j,hgt10j,hgt11j,hgt12j,hgt13j,hgt14j,hgt15j,hgt16j);
HGTpassJaccardDistanceFrame <- data.frame(HGTpassJaccardDistance)
HGTpassJmatrix <- as.matrix.data.frame(HGTpassJaccardDistance)
colnames(HGTpassJmatrix)[1:16] <- c(1:16)
#heatmapJ <- heatmap.plus(HGTpassJmatrix)
heatmapJ <- heatmap(HGTpassJmatrix,symm=TRUE)
plot(heatmapJ) #Jaccard Distance: Measure of DisSimilary (Complement of intersect/union Jaccard Index)


# Need to add better saving
pdf(file="HGTpass-Heatmap-Jaccard-Distance.jpg",height=8,width=8)
plot(heatmapJ)
dev.off()

ggsave(plot(heatmapJ),file="HGTpass-Heatmap-Jaccard-Distance.jpg")

#pdf(file="HGTpass-Heatmap.pdf")
#plot(heatmap)
#dev.off()

##### NEW ANALYSIS 10-30-19 #############

###
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunCA_ugam_bd")
# Ingest vectors for each HGT index and build a new data table/matrix with trees (rows) and a binary 0|1 indicating pass or not for a given HGT constraint (columns)
for (i in HGT){LoadHGTrunVector(i)};
# Build data frame
HGTpass <- cbind(hgt_1,hgt_2,hgt_3,hgt_4,hgt_5,hgt_6,hgt_7,hgt_8,hgt_9,hgt_10,hgt_11,hgt_12,hgt_13,hgt_14,hgt_15,hgt_16);
HGTpass <- data.frame(HGTpass);
colnames(HGTpass)[1:16] <- c(1:16) #numeric column names for each HGT constraint
#write this new table down to disk
write.csv2(HGTpass,"RunCA_ugam_bd_HGTpass.csv");
#build matrix
HGTpassM <-as.matrix(HGTpass)
#Produce heatmap!
heatmap <- heatmap.plus(HGTpassM)
heatmap(HGTpassM)

#to Jaccard Distance (dis-similar): hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt2j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[2],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt3j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[3],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt4j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[4],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt5j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[5],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt6j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[6],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt7j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[7],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt8j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[8],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt9j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[9],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt10j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[10],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt11j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[11],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt12j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[12],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt13j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[13],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt14j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[14],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt15j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[15],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt16j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[16],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );

#Key code to create matrix and heatmap
HGTpassJaccardDistance <- cbind(hgt1j,hgt2j,hgt3j,hgt4j,hgt5j,hgt6j,hgt7j,hgt8j,hgt9j,hgt10j,hgt11j,hgt12j,hgt13j,hgt14j,hgt15j,hgt16j);
HGTpassJaccardDistanceFrame <- data.frame(HGTpassJaccardDistance)
HGTpassJmatrix <- as.matrix.data.frame(HGTpassJaccardDistance)
colnames(HGTpassJmatrix)[1:16] <- c(1:16)
#heatmapJ <- heatmap.plus(HGTpassJmatrix)
heatmapJ <- heatmap(HGTpassJmatrix,symm=TRUE)
plot(heatmapJ) #Jaccard Distance: Measure of DisSimilary (Complement of intersect/union Jaccard Index)

# NEXT ANALYSIS

setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/modeldata/RunCA_cir_bd")
# Ingest vectors for each HGT index and build a new data table/matrix with trees (rows) and a binary 0|1 indicating pass or not for a given HGT constraint (columns)
for (i in HGT){LoadHGTrunVector(i)};
# Build data frame
HGTpass <- cbind(hgt_1,hgt_2,hgt_3,hgt_4,hgt_5,hgt_6,hgt_7,hgt_8,hgt_9,hgt_10,hgt_11,hgt_12,hgt_13,hgt_14,hgt_15,hgt_16);
HGTpass <- data.frame(HGTpass);
colnames(HGTpass)[1:16] <- c(1:16) #numeric column names for each HGT constraint
#write this new table down to disk
write.csv2(HGTpass,"RunCA_cir_bd_HGTpass.csv");
#build matrix
HGTpassM <-as.matrix(HGTpass)
#Produce heatmap!
heatmap <- heatmap.plus(HGTpassM)
heatmap(HGTpassM)

#to Jaccard Distance (dis-similar): hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt1j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[1],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt2j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[2],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt3j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[3],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt4j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[4],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt5j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[5],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt6j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[6],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt7j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[7],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt8j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[8],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt9j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[9],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt10j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[10],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt11j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[11],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt12j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[12],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt13j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[13],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt14j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[14],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt15j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[15],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );
hgt16j <- lapply(seq_along(HGTpass),function(i){ ((1-jaccard(HGTpass[16],HGTpass[i],center=FALSE,px=NULL,py=NULL))*100) } );

#Key code to create matrix and heatmap
HGTpassJaccardDistance <- cbind(hgt1j,hgt2j,hgt3j,hgt4j,hgt5j,hgt6j,hgt7j,hgt8j,hgt9j,hgt10j,hgt11j,hgt12j,hgt13j,hgt14j,hgt15j,hgt16j);
HGTpassJaccardDistanceFrame <- data.frame(HGTpassJaccardDistance)
HGTpassJmatrix <- as.matrix.data.frame(HGTpassJaccardDistance)
colnames(HGTpassJmatrix)[1:16] <- c(1:16)
#heatmapJ <- heatmap.plus(HGTpassJmatrix)
heatmapJ <- heatmap(HGTpassJmatrix,symm=TRUE)
plot(heatmapJ) #Jaccard Distance: Measure of DisSimilary (Complement of intersect/union Jaccard Index)

##### CODE IN DEVELOPMENT #############

#debug
jaccard(hgt_1,hgt_2)
jaccard(HGTpass$hgt_1,HGTpass$hgt_2)
100*(1-jaccard(hgt_1,hgt_2))

##### CODE IN DEVELOPMENT #############

#install.packages('cluster')
library(cluster)
#set.seed()
hgtclusters64<-clara(HGTpass,64,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters32<-clara(HGTpass,32,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters16<-clara(HGTpass,16,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters14<-clara(HGTpass,14,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters12<-clara(HGTpass,12,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters8<-clara(HGTpass,8,metric = c("jaccard"),samples=9900,rngR = TRUE)
hgtclusters4<-clara(HGTpass,4,metric = c("jaccard"),samples=9900,rngR = TRUE)

# example
clara(x, k, metric = c("euclidean", "manhattan", "jaccard"),
      stand = FALSE, samples = 5,
      sampsize = min(n, 40 + 2 * k), trace = 0, medoids.x = TRUE,
      keep.data = medoids.x, rngR = FALSE, pamLike = FALSE, correct.d = TRUE)

# Cluster of interest

hgtcluster <- hgtclusters4

plot(hgtcluster$clustering)
plot(hgtcluster$diss)
plot(hgtcluster$silinfo$widths)


###### CODE IN DEV #########

#Best sample:
#[1]   14  195  841  999 1468 2002 2012 2286 2322 2435 2514 3192 3379 3385 3503 3975 3979 3999 4299 4645 4677
#[22] 4795 4902 5076 5595 5883 5909 6104 6231 6247 6603 6793 6806 7034 7350 7423 7481 7554 7566 7674 7681 8009
#[43] 8018 8127 8448 8607 8627 8631 8659 8698 8808 8820 8867 8868 8899 9072 9162 9365 9369 9506 9633 9675 9690
#[64] 9712 9814 9834 9844 9878

diss <- daisy(hgtclustersmall)
pamhgt <- pam(diss, 16, diss = TRUE)
clusplot(pamhgt,shade = TRUE)

hgt
hgtpam <- pam(HGTpass,8,diss=TRUE,keep.diss = TRUE)
hgtpam$diss
clusplot(hgtpam)

## plotting votes.diss(dissimilarity) in a bivariate plot and
## partitioning into 2 clusters
votes.diss <- daisy(HGTpass)
pamv <- pam(votes.diss, 8, diss = TRUE)
clusplot(pamv, shade = TRUE)
## is the same as
votes.clus <- pamv$clustering
clusplot(votes.diss, votes.clus, diss = TRUE, shade = TRUE)
## Now look at components 3 and 2 instead of 1 and 2:
str(cMDS <- cmdscale(votes.diss, k=3, add=TRUE))
clusplot(pamv, s.x.2d = list(x=cMDS$points[, c(3,2)],
                             labs=rownames(votes.repub), var.dec=NA),
         shade = TRUE, col.p = votes.clus,
         sub="", xlab = "Component 3", ylab = "Component 2")
clusplot(pamv, col.p = votes.clus, labels = 4)# color points and label ellipses
# "simple" cheap ellipses: larger than minimum volume:
# here they are *added* to the previous plot:
clusplot(pamv, span = FALSE, add = TRUE, col.clus = "midnightblue")

