# NetworkFrag_Example.R
# 
# Example of network fragmentation calculations using functions developed by Knight et al (2025)
# for betweeness centrality ranked vertices on a NZ aquaculture network.

# MIT License
# 
# Copyright (c) 2025 Ben Knight, Cawthron Institute
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library(igraph)
library(tibble)
source("robustAnalysis.R") #source the outDegNeigh and robustAnalysis functions

#### Read in data and initialise vars  #########################################################
#read in aquaculture connectivity graph
AgeConn<-readRDS('./data/allconnectionsMSage.rds')  #load connectivity matrices
xy2<-readRDS('./data/allNodelocs.rds') #load vertex locations - in case we need a spatial plot

# targeted removals with out-neighbours
#ageRange= seq(1,30,by=1)  #1:30 # note 2 is for 0-1 days - settlement range time (in days, max = 30 days)
ageRange= c(1,5) # small 1 & 5 day subset of ages to process, can use = 1:30 for all ages
outNeighRng<-0 #out-neighbors to process (set to zero for no out neighbour analysis or 0:4 for up to 4 out neighbours)

#######################################################################################################
################### Process community structure/robustness for given ages and out nbrs ################
######################################################################################################
#### init variables 
nds2plt<-list()
rbLstBtwn<-list()

#loop over all ages
for (a in ageRange){  
    
    print(paste('Loop: ',a))
    age_str=as.character(a) #string index for age
    
    ####  calculate age specific connectivity stats
    conn=AgeConn[,,a]
    colnames(conn)<-1:dim(conn)[1]   #mfishIDs
    conn_trans=conn
    
    # create a graph
    g_transform<-igraph::graph_from_adjacency_matrix(conn_trans, mode = "directed", weighted=T,diag=T) #diag=F - omit self connections?

    #set locations
    igraph::set_vertex_attr(g_transform,"x", value = xy2$x)
    igraph::set_vertex_attr(g_transform,"y", value = xy2$y)

    #calculate number of components and which group each node is part of
    comp=igraph::cluster_infomap(g_transform)
    whichClust<-igraph::membership(comp)
    xy2$memb<-whichClust

    nds2plt[[age_str]]<-data.frame(row.names=as.character(1:dim(xy2)[1]), x=xy2$x, y=xy2$y)
    nds2plt[[age_str]]<-tibble::add_column(nds2plt[[age_str]],btw=0)
    nds2plt[[age_str]]$whichClust[as.numeric(names(whichClust))]<-whichClust
    
    #loop over given range of out neighbours
    for (n in outNeighRng){
      print(paste('Processing out neighbours: ',n))
      outNeighStr<-as.character(n)
     
      # Calc robustness for each age
      rbLstBtwn[[age_str]][[outNeighStr]]<-robustAnalysis(g_transform,numXtra=n)
      nds2plt[[age_str]]$btw<-rbLstBtwn[[age_str]][[outNeighStr]][['nodeImp']]
    }
    print(paste('Age done: ',a))
} #for a  

################################################################################
### Plot robustness metrics for different age and/or out neighbour parameters
##################################################################################
numOutNeighbrs=outNeighRng[1]

a=ageRange[1] #age of network
age_str=as.character(a)
dat1=rbLstBtwn[[age_str]][[1]]

# Example network dismantling plots for chosen "netMetric" metric
netMetric='diam' #set for different metrics chose: "n.comp" or "diam"
nrows1<-dim(dat1$table)[1]
plot(dat1$table$removed.pct[1:(nrows1-1)],
     dat1$table[[netMetric]][1:(nrows1-1)],
     lty=1,type='l',lwd=3,col='gray0',
     xlab='Proportion of nodes removed',ylab = 'Graph diameter',
     main=paste('Fragmentation with node removals, Age:',age_str))
  
#calculate and plot fragmentation index for given age
frag=1-dat1$table$comp.size/dat1$table$nVert #calc fragmentation 
plot(dat1$table$removed.pct[1:(nrows1-1)],frag[1:(nrows1-1)],
     ylim=c(0,1),lty=1,type='l',lwd=3,col='gray0',
     xlab='Proportion of nodes removed',ylab = 'Fragmentation',
     main=paste('Fragmentation with node removals, Age:',age_str)) 


