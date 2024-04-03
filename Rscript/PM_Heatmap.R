#!/usr/bin/env Rscript
#################################################################
# Function: Heatmap
# Call: Rscript PM_heatmap.R -d dist_file [-o outfile]
# R packages used: optparse gplots
# Authors: Yuzhu Chen, Yanhai Gong, Xiaoquan Su
# Updated at July 29, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

## install necessary libraries
p <- c("optparse","pheatmap")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
  make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix file [Required]"),
  make_option(c("-o", "--outfile"), type="character", default='heatmap.pdf', help="Output heatmap [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

## load data
# import dist file
dst_orig <- read.table(file=opts$dist_file, header=TRUE, row.names=1)

## main calc & draw function definition
PM_heatmap <- function(da,clust_method="average") {
  
  dist.data<-as.dist(da)                  # distance matrix to dist structure
  
  #myclust = function(y){                  # function -- hclust with different method
  #  hclust(y,method=clust_method)
  #}
  
  pdf(opts$outfile,width=3+attr(dist.data,'Size')/2,height=3+attr(dist.data,'Size')/2)
  #suppressWarnings(suppressMessages(
  #heatmap.2(as.matrix(dist.data),scale="none",hclustfun=myclust,key=TRUE,
			#key.par=list(cex.axis=2),symkey=FALSE,density.info="none",
			#trace="none",col=colorpanel(100,low="red",high="green"),
            #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
			#lhei=c(1.5, 4, 0.5 ))))
  pheatmap(as.matrix(dist.data),clustering_method = clust_method,border_color = NA,treeheight_row = 200,treeheight_col = 200,angle_col = 45)
  invisible(dev.off())
  
}

## calc
PM_heatmap(da=dst_orig)

