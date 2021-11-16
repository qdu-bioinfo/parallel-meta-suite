#!/usr/bin/env Rscript
#################################################################
# Function: Hcluster
# Call: Rscript PM_hcluster.R -d dist_file [-o outfile -k groupNum]
# R packages used: optparse
# Authors: Yuzhu Chen, Yanhai Gong, Xiaoquan Su
# Updated at July 29, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

## install necessary libraries
p <- c("optparse")
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
  make_option(c("-o", "--outfile"), type="character", default='hcluster.pdf', help="Output file [default %default]"), 
  make_option(c("-c", "--groupNum"), type="integer", default=2, help="Number of groups to rect [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

## load data
# import dist file
dst_orig <- read.table(file=opts$dist_file, header=TRUE, row.names=1)

## main calc & draw function definition
PM_hcluster <- function(da,method="average",k=2) {
  
  dist.data<-as.dist(da)                  # distance matrix to dist structure
  
  hc<-hclust(dist.data,method)            # hclust clustering
  
  pdf(opts$outfile,width=3+attr(dist.data,'Size')/2,height=10)    # need to solve width and height auto adjusting
  par(oma=c(2,2,1,1))
  plot(hc,hang=-1,cex.main=1.6,cex.sub=1.6,cex.lab=1.4,font.sub=1,font.lab=1)
  rect.hclust(hc,k)      # draw
  invisible(dev.off())
  
}

# calc
PM_hcluster(da=dst_orig, k=opts$groupNum)
