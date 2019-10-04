library(ape)
library(coda)
library(mcmcplots)
library(phangorn)
library(plyr)
#### Post sampling the datedist files for HGT constrains
setwd("~/Dropbox (Simons Foundation)/Lily/LM_Trees/")
input = "r1_1.2akinete_sample.datedist" ### insert Model D datedist file here

MyTree = read.tree(input)

conditionalDates = matrix(ncol = length(MyTree[[1]]$node.label),nrow=3000)
conditionalEdges = matrix(ncol = length(MyTree[[1]]$edge.length),nrow=3000)

### Select indices depending on file category
HGTindices = c(24, 3, 27, 83) ## (Sah Transfer donor CF (older), Sah Tranfser recipient CyM (younger),BchH Transfer donor CF (older), BchH Tranfser recipient Ch (younger)

sahDonor = HGTindices[1]
sahRecip = HGTindices[2] 
bchDonor = HGTindices[3]
bchRecip = HGTindices[4]

## post sample the datedist file for tree that meet HGT constraint
n=0
tmp = 0
foo <- vector("list", 2) #create an empty list of 2 elements
class(foo) <- "multiPhylo" #make this list a multiPhylo object
for (i in MyTree){
  tmp = tmp + 1
  dates = as.numeric(i$node.label)
  bchHGT = dates[bchDonor] 
  cyanomel = dates[bchRecip] 
  sahHGT = dates[sahDonor]
  chlorobi = dates[sahRecip]  # crown gsb
  
  if (bchHGT>cyanomel){
    if(sahHGT>chlorobi){
      print(tmp)
      n=n+1
      conditionalDates[n,]=dates
      conditionalEdges[n,]=i$edge.length
      foo[[n]]=i
    }
  }
}


## matrix containing node ages (columns) of trees (rows) that meet the conditional constraint. 
conditionalDates=conditionalDates[complete.cases(conditionalDates),]
conditionalEdges=conditionalEdges[complete.cases(conditionalEdges),]


medPost = apply(conditionalDates,2,FUN = median)
medEdge = apply(conditionalEdges,2,FUN=median)
meanEdge = colMeans(conditionalEdges)
maxPost = apply(conditionalDates,2,FUN = max)
minPost = apply(conditionalDates,2,FUN = min)
ciPost = paste(maxPost,minPost, sep=" - ")

outTable = t(rbind(medPost,ciPost))


newTree = MyTree[[1]]
newTree$node.label=medPost
#newTree$node.label=paste(apply(conditionalDates,2,FUN=max),apply(conditionalDates,2,FUN=min),sep=" - ")

tempTree = MyTree[[1]]
tempTree$node.label = seq(length(tempTree$tip.label)+1,length(tempTree$tip.label)+length(tempTree$tip.label)-1 )
v = rep(0,length(tempTree$edge.length))

edge = tempTree$edge
nodelabels = tempTree$node.label
for (r in 1:nrow(edge)){
  outerNodeIndex = edge[r,1]
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


newTree$edge.length = v


taxonomyKey = read.delim("phototrophy_key_no_bin_key.txt",header=FALSE)
tip = newTree$tip.label
newTree$tip.label=mapvalues(tip, from=taxonomyKey$V3, to=as.character(taxonomyKey$V4))

plot(newTree,show.node.label = FALSE,cex = 0.5,show.tip.label = TRUE)
axisPhylo()


tmp = newTree
tmp$node.label=seq(tmp$node.label)
plot(tmp,show.node.label = TRUE,cex = 0.5,show.tip.label = TRUE)
axisPhylo()
#write.tree(newTree, file="postsampled_chronogram_nobootstrap.tree")
