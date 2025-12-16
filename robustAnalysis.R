# Functions to accompany Knight et al (2025) - Network fragmentation paper.
#
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

# outDegNeigh(graphIn,numXtra =0 ,val)
#
# Function to find, up to numextra, highest importance out-degree neighbors for the most
# important vertex in a graph.
#
# Inputs:
# graphIn = igraph object with named vertices and a connected directional graph
# numXtra = maximum number of extra out degree neighbours
# val = importance  values for each vertex in the graph (e.g. betweeness centrality values)
#
# Output: names of neighbour vertices
# Example:
#    g<-make_tree(10) #make a random tree graph
#    V(g)$name<-as.character(1:length(V(g))) #set names based on orignal vertex order
#    E(g)$weight<-runif(length(E(g))) #give g random edge weights (connection probabilities)
#    out<-outDegNeigh(g,2,val)
outDegNeigh<-function(g,numXtra,val){
  #check for valid input graph
  stopifnot(igraph::is_igraph(g));
  
  # rank graph vertices by their importance
  impOrder<-igraph::V(g)$name[order(val, decreasing = TRUE)];
  
  # find most important 1st degree connected out neighbour vertices
  connImp<-igraph::ego(graph = g,order = 1,nodes = impOrder[1],mode = 'out')[[1]];
  
  # order by highest importance node
  btwCentrOutOrdr<-order(val[connImp], decreasing = TRUE);
  
  # output any extra out neighbour vertices 
  return(impOrder[btwCentrOutOrdr[1:(numXtra+1)]]);
}

# out <- robustAnalysis(g, numXtra=0, impmethod = "btwn.cent", val=NA)
# Function to estimate variation in robustness metrics for a graph (g) from 
# structured node removals of high value  to low value (val) vertices
#
# Inputs:
# g = igraph object of connection probabilities, specified as edge "weight"
# numXtra = number of extra out-degree neighbors to remove with each targeted removal (0=none)
# impmethod = method to determine importance of nodes for removal, 
#                  For impmethod specify one of:
#                  "btwn.cent" (betweenness centrality), 
#                  "degree" (degree centrality), 
#                 "closenness" (closeness centrality), or 
#                  "manual" (use val supplied), if impmethod = "manual", val variable must be supplied
# val [Optional, except for when impmethod is "manual"] 
#            val is used to determine the importance order of the vertices (high value = high importance) 
#            Can be calculated from metrics such as: betweenness centrality
#            If "val" is not provided, impmethod will be used to determine values
# Output: Returns a list output with node importance values (nodeImp), and 
#          a table containing columns with:
#            comp.size = number of vertices in the largest community/cluster
#            comp.meansize = mean number of vertices in all communities/clusters
#            comp.n = number of communities/clusters
#            diam = diameter of the graph [mean distance log10(1/connection probabilities)]
#            comp.pct = fraction of graph vertices remaining
#            removed.pct = fraction of graph vertices removed
#            nVert = graph size - number of vertices remaining
# Examples: 
#           g<-make_tree(10) #make a random tree graph
#           V(g)$name<-as.character(1:length(V(g))) #set names based on orignal vertex order
#           E(g)$weight<-runif(length(E(g))) #give g random edge weights (connection probabilities)
#           out<-robustAnalysis(g, numXtra=0, impmethod = "btwn.cent") #run robustness analysis 
#           plot(g) #view the graph?
#          
robustAnalysis<-function (g, numXtra=0, impmethod = "btwn.cent", val=NA){ 
  
  #check for valid input graph
  stopifnot(igraph::is_igraph(g));
  impmethod=tolower(impmethod); #make lower case
  
  #create temporary graph "g_dist" with distance based metrics from connectivity probabilities for centrality calculations
  g_dist=g
  badval=30; # large exponent value to set up small non-zero distance
  dist<-log10(1/(igraph::E(g)$weight+10^-badval))  #note add small value to stop log zero issues
  igraph::E(g_dist)$weight<-dist
  badE=igraph::E(g_dist)$weight==badval 
  g_dist<-igraph::delete_edges(g_dist, which(badE))  #now delete bad edge weights (no connections)
  
  cmp<-igraph::cluster_infomap(g); #calculate infomap community clusters
  
  #get infomap community/cluster metrics
  csize<-igraph::sizes(cmp);
  orig_max <- max(csize);
  orig_mean<- mean(csize);
  orig_comp<- length(csize) 
  orig_d<-igraph::diameter(g);
  
  #initialise outputs
  n <- igraph::vcount(g)  
  cumremoved <- seq.int(0, 1, length.out = n + 1L);
  max.comp.removed <- rep.int(orig_max, n);
  mean.comp.removed <- rep.int(orig_mean, n);
  n.comp <- rep.int(orig_comp, n);
  diam<-rep.int(0, n);
  out<-list()
  
  if (impmethod != "manual") {
    if (impmethod == "btwn.cent"){
      val <- igraph::centr_betw(g_dist)$res;
    } else if(impmethod == "degree"){
      val <- igraph::degree(g_dist)
    } else if(impmethod == "closenness"){
      val <- igraph::closenness(g_dist)
    }else if(impmethod == "manual"){ 
      if (any(is.na(val))){stop("Vertex importance values ('val') required for impmethod = 'manual'")} 
    }else{
      stop(paste("Unrecognised value for 'impmethod':", impmethod))  
    }
    ord <- igraph::V(g)$name[order(val, decreasing = T)];
  } #if impmethod
  
  nVert<-rep.int(length(ord), n); #init vertice outputs
  out[['nodeImp']]<-val  #save original importance metrics
  
  j=1; #init loop counter
  cumDel=0; #init removal counter
  
  while (length(g) > numXtra) {
     #find highest ranking node (and highest out neighbours if numXtra >0)
    delVert<- outDegNeigh(g,numXtra,val); 
    keep<-!is.na(delVert) #check for NAs
    
    #note val inds are different to vertex names
    valDel<-which(!is.na(match(igraph::V(g)$name,delVert[keep])))
    
    ##reduce graph size
    g <- igraph::delete_vertices(g, delVert[keep]);
    g_dist<-igraph::delete_vertices(g_dist, delVert[keep]);
    val<- val[-valDel] #reduce number of vals
    
    #update network metrics for reduced graph
    if(length(g) > numXtra){ 
      cmp<-igraph::cluster_infomap(g);
      csize<-igraph::sizes(cmp)
      mean.comp.removed[j+1L] <- mean(csize,na.rm = T);
      max.comp.removed[j+1L] <- max(csize);
      nVert[j+1L] <- length(igraph::V(g));
      n.comp[j+1L] <- length(csize);
      diam[j+1L] <- igraph::diameter(g_dist);
      cumremoved[j+1L] <-cumDel
      j=j+1
      cumDel<-cumDel+sum(keep);
    } else{
      break
    }
  }
  
  # compile output metrics for analysis and plotting
  nVert<-c(nVert[2:j],0)
  cumremoved<-c(cumremoved[2:j],0)
  removed.pct <- cumremoved/cumDel;
  max.comp.removed <- c(max.comp.removed[2:j],0)
  mean.comp.removed <- c(mean.comp.removed[2:j], 0)
  n.comp<- c(n.comp[2:j], 0)
  diam<- c(diam[2:j], 0)
  
  # calculate percet removed as a fraction of original
  comp.pct <- max.comp.removed/orig_max
  
  #output metric changes as a summarized table 
  out[['table']] <- data.table::data.table(comp.size = max.comp.removed, 
                                           comp.meansize = mean.comp.removed, comp.n=n.comp, diam=diam,
                                           comp.pct = comp.pct, removed.pct = removed.pct,
                                           nVert=nVert)
  return(out)
}  #robustAnalysis
