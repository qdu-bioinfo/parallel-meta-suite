#!/usr/bin/env Rscript
#################################################################
# Function:  Co-occurence Network analysis
# Call: Rscript PM_network.R -i dist_file -o outfile -t threshold
# R packages used: optparse
# Last update: 2016-09-29, Xiaoquan Su, Honglei Wang
#            : The color and size of each nodes could change with
#            : its connective numbers. 
#################################################################

## install necessary libraries
p <- c("optparse","igraph" )
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')
##
map <- function(x, range = c(0,1), from.range=NA) {
    if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)
    
    ## check if all values are the same
    if(!diff(from.range)) return(
	    matrix(mean(range), ncol=ncol(x), nrow=nrow(x), 
		    dimnames = dimnames(x)))
    
    ## map to [0,1]
    x <- (x-from.range[1])
    x <- x/diff(from.range)
    ## handle single values
    if(diff(from.range) == 0) x <- 0 
    
    ## map from [0,1] to [range]
    if (range[1]>range[2]) x <- 1-x
    x <- x*(abs(diff(range))) + min(range)
    
    x[x<min(range) | x>max(range)] <- NA
    
    x
}


## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
option_list <- list(
  make_option(c("-i", "--dist_file"), type="character", help="Input distance matrix [Required]"),
  make_option(c("-o", "--outfile"), type="character", default='network.pdf', help="Output Network [default %default]"),
  make_option(c("-p", "--positive_edges"), type="logical", default=T,help="If enable the positive edges [Optional, default %default]"),
  make_option(c("-n", "--negative_edges"), type="logical", default=T,help="If enable the negative edges [Optional, default %default]"),
  make_option(c("-t", "--threshold"), type="double", default=0.7, help="Edge threshold [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

Corr <- read.table(opts$dist_file, header=T, row.names=1)
Threshold <- opts$threshold
Pos_Edge <- opts$positive_edges
Neg_Edge <- opts$negative_edges

#init graph
if (Pos_Edge & Neg_Edge){
   g <- graph.adjacency(abs(Corr) > Threshold & abs(Corr)<1,mode="undirected")
}else if (Pos_Edge){
     g <- graph.adjacency(Corr > Threshold & Corr < 1,mode="undirected")
}else if (Neg_Edge){
     g <- graph.adjacency(Corr < -1 * Threshold & Corr > -1 ,mode="undirected")
}else stop('Please check the pos/neg edges parameters')

E(g)$color <- 'grey'
V(g)$size <- 10
V(g)$color <- 'cornflowerblue'
E(g)$label <- 1

Pos_Edge_N <- 0
Neg_Edge_N <- 0

for(i in 1:(length(Corr)-1)){
    for(j in i:length(Corr)){
		if ((Neg_Edge) & (Corr[i,j] < -1 * Threshold) & (Corr[i,j] > -1)){
			E(g, path=c(colnames(Corr)[i],colnames(Corr)[j]))$color <- 'red'
			E(g, path=c(colnames(Corr)[i],colnames(Corr)[j]))$weight <- round(abs(Corr[i,j]),3)
			Neg_Edge_N <- Neg_Edge_N + 1
		}
		if ((Pos_Edge) & (Corr[i,j] > Threshold) & (Corr[i,j] < 1)){
			E(g, path=c(colnames(Corr)[i],colnames(Corr)[j]))$color <- 'green'
			E(g, path=c(colnames(Corr)[i],colnames(Corr)[j]))$weight <- round(abs(Corr[i,j]),3)
			Pos_Edge_N <- Pos_Edge_N + 1
		}	
	}
}


#****PDF name****
pdf(opts$outfile,width=15,height=15)
str = paste("Threshold: ",Threshold,"  Nodes: ",length(V(g)),"  Islands: ",clusters(g)$no, "\n# of Pos edge: ",Pos_Edge_N, "  # of Neg edge: ",Neg_Edge_N, "  Edge length: ",average.path.length(g),"\nDensity: ",graph.density(g),"  Diameter: ",diameter(g),"  Radius: ",radius(g),"  Centralization:",centralization.closeness(g)$centralization)
plot(g, layout=layout.fruchterman.reingold,main=str,edge.label=E(g)$weight, vertex.size=map(degree(g),c(1,10)), vertex.color=map(degree(g),c(1,10)))
legend("top",c("pos","neg"),lty=c(1,1),col=c("green","red"))
invisible(dev.off())
