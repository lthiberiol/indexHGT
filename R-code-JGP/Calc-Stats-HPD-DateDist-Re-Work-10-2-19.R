################################################################################################
### Massachusetts Institute of Technology                                                      #
### Department of Earth, Atmospheric & Planetary Sciences, Fournier Lab                        #
### Histogram-Node-script                                                                      #
### Purpose: Statistics for Posterior distribution given datedist | node                       #
### Authors: Jack G. Payette, M.S. (https://orcid.org/0000-0001-6479-4557)                     #
### Authors2: Professor Gregory P. Fournier, Ph.D.                                             #
### Contribution: John K. Kruschke: http://dx.doi.org/10.1016/B978-0-12-405888-0.00025-8       #
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

################################################################################################
### BEGIN SCRIPT EXECUTION #####################################################################
################################################################################################

### Set Working Directory from MIT Dropbox ###
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/")
#setwd("Users/jackpayette/Documents/GitHub/R-code-mit") #Home-Mac WorkDir
#setwd("/Users/Jack Payette/Dropbox (MIT)/R-code-mit/R-code-mit/") #Home-PC WorkDir
print(getwd())

################################################################################################
### DEFINE FUNCTIONS ###########################################################################
################################################################################################
# Work on putting these in separate source("NAMEHERE.R")

### LoadDateDist Function: Example underneath.
LoadDateDist <-function(DateDistFilePath,UniqueModelCode){
  model <- read.tree(DateDistFilePath);
  class(model) <- "multiPhylo";
  assign(paste0(UniqueModelCode),model,envir = .GlobalEnv);
  assign(paste0(UniqueModelCode,"datedist",sep=""),DateDistFilePath,envir = .GlobalEnv);
} 
# LoadDateDist("modeldata/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB") 

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
}
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
  returnValue(NodeAge)}
# PrintRangerNodeAge(MyTreeLabels,MyTreeLabelsRanger,BB,1,'n61') 

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

### Statistics Function: Highest Posterior Density Interval | DateDistNode Vector Object
HPDIstats = function(sampleVec,credMass=0.95){
  # Computes highest density interval from a sample of representative values,
  # estimated as shortest credible interval.
  # Arguments:
  # sampleVec
  # is a vector of representative values from a probability distribution.
  # credMass
  # is a scalar between 0 and 1, indicating the mass within the credible
  # interval that is to be estimated.
  # Value:
  # HDIlim is a vector containing the limits of the HDI sortedPts = sort(sampleVec)
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling(credMass * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0, nCIs)
  for(i in 1:nCIs) {
    ciWidth[i]= sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  
  PosteriorStdDev<-sd(sampleVec)
  PosteriorMean<-mean(sampleVec)
  PosteriorMedian<-median(sampleVec)
  PosteriorSampleN<-length(sampleVec)
  
  marginError <- qnorm(((1-credMass)/2)+credMass)*((sd(sampleVec))/sqrt(length(sampleVec)))
  LowerConfidenceInterval = PosteriorMean - marginError
  UpperConfidenceInterval = PosteriorMean + marginError
  
  HDI_stats <- cbind(HDImin,HDImax,PosteriorMedian,PosteriorSampleN,PosteriorMean,PosteriorStdDev,marginError,LowerConfidenceInterval,UpperConfidenceInterval)
  HDI_stats <- data.frame(HDI_stats)
  return(HDI_stats)
}
HPDIstats(BB177, credMass=0.95)

################################################################################################
### DATA             ###########################################################################
################################################################################################

#### Load DateDist | Node
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",329,"chlorobiales crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",293,"chloroflexia crown")
## Label DataDist x Node
Cyanobacteria<-BB177
Chlorobiales<-BB329
Chloroflexia<-BB293

### Set Working Directory from MIT Dropbox ###
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/DraftSpace/")
### Write Sampled data down into CSV
write.csv(Cyanobacteria,"BBcyanobacteriaHistogram.csv")
write.csv(Chlorobiales,"BBchlorobialesHistogram.csv")
write.csv(Chloroflexia,"BBchloroflexiaHistogram.csv")
### Basic density plots
plot(density(Cyanobacteria),xlim=c(5000,0))
plot(density(Chlorobiales),xlim=c(5000,0))
plot(density(Chloroflexia),xlim=c(5000,0))


PhyloBayesModel <- cbind(Cyanobacteria,Chlorobiales,Chloroflexia); PhyloBayesModel <- data.frame(PhyloBayesModel);
PhyloBayesModelMelt <- melt(PhyloBayesModel, varnames=c(Cyanobacteria,Chlorobiales,Chloroflexia))

#PhyloBayesModel <- cbind(Cyanobacteria); PhyloBayesModel <- data.frame(PhyloBayesModel);
#PhyloBayesModelMelt <- melt(PhyloBayesModelMelt, varnames=c(Cyanobacteria))

# Write melted data file down
write.csv(PhyloBayesModelMelt,"Data-model-BB-Nodes-Age-177-329-293.csv")

# Basic plot for melted data
ggplot(PhyloBayesModelMelt, aes(value, col = variable)) + scale_x_reverse(limits=c(5000,0)) + geom_freqpoly() 
# Melted Data Plot
gplot<- ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank());
plot(gplot)
#Enhanced statistics 
g <- print(gplot)
head(g)
write.csv(data.frame(g$data),"Density-plot-model.csv")
#Print function
print(get("compute_group", ggplot2::StatDensity))


#######################3

##### FRESH LOAD ALL MODELS
### Set Working Directory from MIT Dropbox ###
setwd("/Users/payette/Dropbox (MIT)/R-code-mit/R-code-mit/")

#model	DateDistFile
#### Load DateDist | Node
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBA_ugam_bd_7_20_sample.datedist","BA",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBC_ugam_bd_7_20_sample.datedist","BC",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBD_ugam_bd_7_20_sample.datedist","BD",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBE_ugam_bd_7_20_sample.datedist","BE",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBF_ugam_bd_7_20_sample.datedist","BF",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","modeldata/Cyano_modelBG_ugam_bd_7_20_sample.datedist","BG",177,"cyanobacteria crown")

PhyloBayesModel <- cbind(BA177,BB177,BC177,BD177,BE177,BF177,BG177); PhyloBayesModel <- data.frame(PhyloBayesModel);
PhyloBayesModelMelt <- melt(PhyloBayesModel, varnames=c(BA177,BB177,BC177,BD177,BE177,BF177,BG177))

gplot<- ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,2000)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank());
plot(gplot)
#Enhanced statistics 
g <- print(gplot)
head(g)
write.csv(data.frame(g$data),"Density-plot-model.csv")
#Print function
print(get("compute_group", ggplot2::StatDensity))

#################################### NEW PLOT BELOW

ggsave("Density-plot-7-models-Nodes-Age-177.PDF",
       ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
       ) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                        bw = "nrd0", adjust = 1,
                        kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                        show.legend = NA, position="identity"
       )+scale_x_reverse(limits=c(4000,2000)
       )+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
              caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
       )+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())
       ,units = "in",scale = 1,height=6,width=8)

gplot<- ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,2000)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank());
plot(gplot)
#Enhanced statistics 
g <- print(gplot)
head(g)
write.csv(data.frame(g$data),"Density-plot-model.csv")
#Print function
print(get("compute_group", ggplot2::StatDensity))

#################################### NEW PLOT ABOVE


####### IN DEVVVVVVVVVV

BA<-HPDIstats(BA177, credMass=0.95)
BA$LowerConfidenceInterval
BA$UpperConfidenceInterval

BA <- HPDIstats(BA177, credMass=0.95)
BB <- HPDIstats(BB177, credMass=0.95)
BC <- HPDIstats(BC177, credMass=0.95)
BD <- HPDIstats(BD177, credMass=0.95)
BE <- HPDIstats(BE177, credMass=0.95)
BF <- HPDIstats(BF177, credMass=0.95)
BG <- HPDIstats(BG177, credMass=0.95)

BA$HDImin
BA$HDImax
BA177
PhyloBayesModel <- cbind(BA177); PhyloBayesModel <- data.frame(PhyloBayesModel);
PhyloBayesModelMelt <- melt(PhyloBayesModel, varnames=c(BA177))


BA<- data.frame(BA177,c(1:100))

############ IN DEVVVVVVVVVVV
ggplot(PhyloBayesModelMelt,aes(value, shape=variable)
) +scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank()
)+geom_linerange(aes(y=PhyloBayesModelMelt$value,ymin=BA$LowerConfidenceInterval, ymax=BA$UpperConfidenceInterval), width=50)

##
ggplot(PhyloBayesModelMelt,aes(value, shape=variable)
) + geom_bar(aes(value),alpha = 0.1,inherit.aes = FALSE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank()
)+geom_errorbar(aes(x=PhyloBayesModelMelt$value,ymin=BA$HDImin, ymax=BA$HDImax), width=40)

##
ggplot(PhyloBayesModelMelt,aes(value, color=variable)
) +scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from 7 different Phylobayes model runs (Nodes#177).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank()
)+geom_errorbar(aes(x=PhyloBayesModelMelt$value,ymin=BA$HDImin, ymax=BA$HDImax), width=40)+ stat_density(alpha = 0.1,
               bw = "nrd0", adjust = 1,
               kernel = "gaussian",trim = FALSE, na.rm = FALSE,
               show.legend = NA, position="identity"
)
par(new=TRUE)
ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + geom_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())
#
ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + geom_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())+geom_errorbar(aes(ymin=BA$HDImin, ymax=BA$HDImax), width=40)







### IN DEV FUNCTIONS ###########################################################################



##############
##### from Plot script
PlotHistModelNode<-function(UniqueModelCode,NodeNumberToFilter,ModelNode){
  #Can't get this to work yet: assign("ModelNode",as.name(paste0(UniqueModelCode,NodeNumberToFilter)),envir=.GlobalEnv)
  ModelTrees<-length(ModelNode)
  pdf(paste0(paste0(paste0(paste0("Histogram-model-",UniqueModelCode),"-Nodes-Age-"),NodeNumberToFilter),".pdf"), width=10, height=6) 
  hist(ModelNode,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",
       ylab=paste0(paste0("Percent of Total (Frequency; N=",ModelTrees),")"),main=(paste0(paste0(paste0("Histogram of Age Estimates for Node#",NodeNumberToFilter)," for PhyloBayes model:\n"),UniqueModelCode)));
  dev.off();
}
 PlotHistModelNode("BB",177,BB177)

###################################
### Initial dev above complete
### Debugging below:
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",177,"cyanobacteria crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",329,"chlorobiales crown")
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",293,"chloroflexia crown")
########
LoadDateDistNodeComments("data/Cyano_modelBC_ugam_bd_7_20_sample.labels","data/Cyano_modelBB_ugam_bd_7_20_sample.datedist","BB",177,"cyanobacteria crown")
hist(BB177,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",ylab=paste0(paste0("Percent of Total (Frequency; N=",NumberOfTreesInBB),")"),
     main=paste0(paste0(paste0("Histogram of Age Estimates for Node#",NodeNumberToFilter)," for PhyloBayes model:\n"),BBDateDist))
##########################

######

pdf("data/GSA/Histogram-model-BB-Nodes-Age-177-329-293.pdf", width=10, height=6) 
par()
histogram(Cyanobacteria,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",ylab="Percent of Total (Frequency; N=100)",main="Histogram of Age Estimates for Cyanobacteria Crown (Node#177) for PhyloBayes model BB.")
histogram(Chlorobiales,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",ylab="Percent of Total (Frequency; N=100)",main="Histogram of Age Estimates for Chlorobiales Crown (Node#329) for PhyloBayes model BB.")
histogram(Chloroflexia,xlim=c(5000,0),ylim=c(0,50),xlab="Age in Ma (Millions of Years Ago)",ylab="Percent of Total (Frequency; N=100)",main="Histogram of Age Estimates for Chloroflexia Crown (Node#293) for PhyloBayes model BB.")
dev.off()

#node 177 (cyano crown) BBcyanobacteria
#node 329 (chlorobiales crown) BBchlorobiales
#node 293 (chloroflexia crown) BBchloroflexia

write.csv(Cyanobacteria,"data/GSA/BBcyanobacteriaHistogram.csv")
write.csv(Chlorobiales,"data/GSA/BBchlorobialesHistogram.csv")
write.csv(Chloroflexia,"data/GSA/BBchloroflexiaHistogram.csv")

PhyloBayesModel <- cbind(Cyanobacteria,Chlorobiales,Chloroflexia)
PhyloBayesModel <- data.frame(PhyloBayesModel)
PhyloBayesModelMelt <- melt(PhyloBayesModel, varnames=c(Cyanobacteria,Chlorobiales,Chloroflexia))
PhyloBayesModelMelt <- melt(PhyloBayesModelMelt, varnames=c(Cyanobacteria))

# Write melted data file down
write.csv(PhyloBayesModelMelt,"data/GSA/Data-model-BB-Nodes-Age-177-329-293.csv")

# Basic plot for melted data
ggplot(PhyloBayesModelMelt, aes(value, col = variable)) + scale_x_reverse(limits=c(5000,0)) + geom_freqpoly() 
# Basic density plots
plot(density(Cyanobacteria),xlim=c(5000,0))
plot(density(Chlorobiales),xlim=c(5000,0))
plot(density(Chloroflexia),xlim=c(5000,0))

ggsave("data/GSA/Density-plot-model-BB-Nodes-Age-177-329-293.PDF",
       ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
       ) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                        bw = "nrd0", adjust = 1,
                        kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                        show.legend = NA, position="identity"
       )+scale_x_reverse(limits=c(4000,0)
       )+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
              caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
       )+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())
       ,units = "in",scale = 1,height=6,width=8)
#New adjustment to use 6x8 inch for PDF

# plot
ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())

#store stats data from plot fn

gplot<-               ggplot(PhyloBayesModelMelt,aes(value, col=variable,fill=variable)
) + stat_density(alpha = 0.1,inherit.aes = TRUE,
                 bw = "nrd0", adjust = 1,
                 kernel = "gaussian",trim = FALSE, na.rm = FALSE,
                 show.legend = NA, position="identity"
)+scale_x_reverse(limits=c(4000,0)
)+labs(title="Density Plot for Age Estimates of Crown Cyanobacteria, Chlorobiales & Chloroflexia.",
       caption = paste0("\n N.B. Fixed X-axis reversed in Ma. R Function: stat_density, bw = nrd0, kernel = gaussian.\n Data Source: 100 Chronogram Trees from Phylobayes model run BB (Nodes#177,329,293).\n Fournier Lab, MIT.")
)+xlab("Age in Ma (Millions of Years Ago)")+theme(legend.title = element_blank())


#Enhanced statistics 
g <- print(gplot)
head(g)
write.csv(data.frame(g$data),"data/GSA/Density-plot-model-BB-Nodes-Age-177-329-293.csv")
#Print function
print(get("compute_group", ggplot2::StatDensity))
write.csv(get("compute_group", ggplot2::StatDensity),"data/GSA/Density-formula.csv")

##############################################
###                      #####################
##############################################

################################################################################################
### DEFINE FUNCTIONS ###########################################################################
################################################################################################

### FROM HGT CONSTRAINTS SCRIPT

### Filter/Post sample a given datedist for ONLY Trees that pass HGT constraints 
#Inputs: - multiPhylo datedist file = MyTree ***
#        - HGT constraint for internal node # ordered index
#        - a,d are passed in from HGTnode_input file which maps to tree labels file
# *** Ancestor node OLDER than Descendant node; a > d ***
#        a row in HGT index file = run name
#Outputs:a "PassedTrees.datedist" with prefix of run name.
#TO-DO FIX ERROR ON <2 TREES IN LOOP!!!
#TO-DO allow a project name to be passed into this function and written to output files
### FIXED FUNCTION 9-30-19
SampleTrees <- function(tree0,a,d,nameofrun,print=T){
  MyTree <- tree0
  n=0
  tmp = 0
  nPassedTrees =0
  PassedTrees <- vector("list", 2) #create an empty list of 2 elements to store PassedTrees
  class(PassedTrees) <- "multiPhylo" #make this list a multiPhylo object
  write.tree(MyTree,"AllInitialTrees.datedist")
  AllTrees <- read.tree(file="AllInitialTrees.datedist")
  class(AllTrees) <- "multiPhylo"
  TreePassedConstraint<-as.vector(c(seq(AllTrees)))
  #Dates
  conditionalDates = matrix(ncol = length(AllTrees[[1]]$node.label),nrow=3000)  #(optional) a vector of mode character giving the names of the nodes.
  conditionalEdges = matrix(ncol = length(AllTrees[[1]]$edge.length),nrow=3000)  #(optional) a numeric vector giving the lengths of the branches given by edge.
  # For loop to evaluate one HGT
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
    }
  }
  #Write PassedTrees given A>D in Newick format from ape in R... NEED TO CHECK APPROPRIATE FILE FORMATTTING
  if ( (nPassedTrees<2)==FALSE){
    filename1 <- "PassedTrees.datedist"
    filename1 = paste0(nameofrun,filename1,sep="")
    write.tree(PassedTrees,file=filename1)
    filename2 <- "PassedTreesIndex.txt"
    filename2 = paste0(nameofrun,filename2,sep="")
    write.csv(TreePassedConstraint,file=filename2)
    return(PassedTrees)} else {return("< 2 Trees satisfied by HGT Constraint Error. Script Stopped.")}
};


### Loop over all HGT constraints and for each, save PassedTrees in a new datedist file with prefix# as row of HGT constraint.
#TO-DO: Have better failure/error handling for returning 0 or 1 trees!!!
#NOTE TO-DO: Include the Node naming checks within this function!
#LOOP FOR each i in HGTnode
#Inputs: a multiPhylo datedist file = MyTree = tree0
#        a HGT index  
#        a ProjectName = nameofrun
#Outputs: to console, table/Text
#Outputs: DATA from sub function called SampleTrees() outputs to datedist file PassedTrees.txt
FilterTreeHGT <- function(tree2,hgts,nameofrun,print=T)
{
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
    SampleTrees(MyTree,a,d,i) #Given a set of Trees and constraints, filter and write results to file
    print(paste(i,c("HGT constraint run complete!"),sep=" "))
    PassedHGT[[n]]=  ReadPassedTrees(i) #Not the cleanest results function/table but can adjust.
  }
  write.table(PassedHGT,paste(nameofrun,"PassedHGTResults.txt"))
  return(PassedHGT) #Need to fix since this returns blank!
};

### Check & Fail for mis match in number of internal nodes between trees in original datedist and from Reference Label Tree
#NOTE TO-DO: Have a way of checking order not just node number
#Inputs: a multiPhylo datedist file = MyTree = tree1
#        a Reference Label Tree from file = MyTreeLabels = tree2
#Output: Returns NULL if Passed, Returns Failed Tree i from datedist if Assumption not met
CheckFailInternalNodeN <- function(tree1,tree2,print=T) #Given tree 2 is reference label
{
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
};

# Read a given datedist file from an HGT filter/sample run and compare to original datedist and output
# Input: a run number corresponding to HGT index
# Outputs: text file
ReadPassedTrees <- function(nameofrun,print=T)
{
  AllTrees <- read.tree(file="All Initial Trees.datedist")
  class(AllTrees) <- "multiPhylo"
  #Create new function to load back trees from output to test
  filename1 <- "PassedTrees.datedist"
  filename1 = paste(nameofrun,filename1)
  PassedTrees1 <- read.tree(file=filename1)
  checkoutput <-paste(paste("Results:",paste(length(PassedTrees1),length(AllTrees),sep=" of ")),"Trees passed HGT constraints specified from datedist & file inputs.")
  return(checkoutput)
}

# Check a given datedist file from an HGT filter/sample run
# Inputs: a run number corresponding to HGT index
#         a Refernce Label Tree = MyTreeLabels = tree2
CheckPassedTrees <- function(nameofrun,tree2,print=T)
{
  filename1 <- "PassedTrees.datedist"
  filename1 = paste(nameofrun,filename1)
  PassedTrees <- read.tree(file=filename1)
  class(PassedTrees) <- "multiPhylo"
  Cresult <- CheckFailInternalNodeN(PassedTrees,tree2)
  checkoutput <- is.null(Cresult)
  return(checkoutput)
}


# OLD Functions
### IN DEV FUNCTIONS ###########################################################################
#Utility Fn # Requires separate MyTreeLabels call
PrintNodeAge <- function(MyTreeInput,TreeNumber,StandardNodeNumber){
  Nnodes <- MyTreeInput$Nnode
  internalNodeLabels = data.frame(N = c(seq(1:Nnodes)),Node=c((MyTreeInput$node.label)))
  NodeAge<-as.numeric(MyTreeInput[[TreeNumber]]$node.label[with(internalNodeLabels,N[Node==(StandardNodeNumber)])])
  returnValue(NodeAge)}

# From original script to plot the median node age on a single tree from a given datedist file
# Writes to a table and PDF
PlotPassedTrees <- function(datedistfile)
{
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

# From original script to plot the median node age on a single tree from a given datedist file
# Writes to a table and PDF
PlotPassedTrees <- function(datedistfile)
{
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

##############################################
### END SCRIPT EXECUTION #####################
##############################################
  
