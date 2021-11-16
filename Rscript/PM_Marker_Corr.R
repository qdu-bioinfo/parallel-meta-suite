#!/usr/bin/env Rscript
#################################################################
# Function: Correlation analysis between taxa and each of numerical metadata
# Call: Rscript PM_Marker_Corr.R -i abund_file -m mapfile(metadata) -o outfile
# R packages used: psych and optparse
# Authors: Yuzhu Chen, Shi Huang, Zheng Sun, Xiaoquan Su
# Updated at Aug. 20, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################
#-----------install necessary libraries--------------------------------
p <- c("optparse","psych")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
  suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))
#---------------------------------------------------------------------
## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--table_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),
  make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
  make_option(c("-o", "--out_dir"), type="character", default='Corr_features', help="Output file path [default %default]"),
  make_option(c("-p", "--prefix"), type="character", default='Out', help="Output file prefix [default %default]"),
  make_option(c("-t", "--p_cutoff"), type="double", default=0.1, help="The cutoff of adjusted P values [default %default]"),
  make_option(c("-r", "--r_cutoff"), type="double",default=0.4, help="The cutoff of correlation coefficients [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$table_file)) stop('Please input a feature table (*.Abd)')
if(is.null(opts$meta_data)) stop('Please input a meta data file')

# load data
filename <- opts$table_file
metadata.filename <- opts$meta_data
outpath <- opts$out_dir
prefix <- opts$prefix
p.adj.cutoff <- opts$p_cutoff
r.cutoff <- opts$r_cutoff

#---------------------------------------------------------------------
dir.create(outpath)
g<-read.table(filename,header=T,row.names=1,sep="\t"); g<-g[order(rownames(g)),]
gmat<-data.matrix(g)
#---------------------------------------------------------------------
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1,as.is=FALSE); allmetadata<-allmetadata[order(rownames(allmetadata)),]
q_metadata<-data.frame(allmetadata[sapply(allmetadata,class)!="factor"][order(rownames(allmetadata)),])
f_metadata<-data.frame(allmetadata[sapply(allmetadata,class)=="factor"][order(rownames(allmetadata)),])
names(q_metadata)<-q_metadata.name<-names(which(sapply(allmetadata,class)!="factor"))
q_metadata.num<-length(which(sapply(allmetadata,class)!="factor"))
#---------------------------------------------------------------------
#if(q_metadata.num==0) { cat('Warning: No numerical data contained in the meta data file, PM_Marker_Corr exit\n')}else{
if(q_metadata.num!=0) {
#---------------------------------------------------------------------
#cat("The number of variables: ", ncol(g) ,"\n")
#cat("The number of metadata: ", ncol(allmetadata), "\n")

corr_m<-corr.test(g,q_metadata,method="spearman")
r<-corr_m$r
p<-corr_m$p
#sink(paste(outpath,"/", prefix, ".Corr.all.r.xls",sep=""))
#cat("\t"); write.table(r,quote=FALSE,sep="\t"); sink()
#sink(paste(outpath,"/",prefix, ".Corr.all.pvalues.xls",sep=""))
#cat("\t"); write.table(p,quote=FALSE,sep="\t"); sink()

#--------------------------------------------------
# ScatterPlot
#--------------------------------------------------
for(a in 1:q_metadata.num){
   #------------------Plot for each of metadata
   #cat("Bacterial correlation with", q_metadata.name[a],"\n")
    if(q_metadata.num>1){
        out<-data.frame(R=r[,a], Adj_P=p[,a])
        corr_ord<-with(out,order(-abs(R)))
        out<-out[corr_ord,]
        mat<-gmat[,corr_ord]
        sink(paste(outpath,"/", prefix, ".", q_metadata.name[a],".Corr.xls",sep="")); cat("\t"); write.table(out,quote=FALSE,sep="\t");sink()
     }else{     
        out<-data.frame(R=r, Adj_P=p)
        corr_ord<-with(out,order(-abs(R)))
        out<-out[corr_ord,]
        mat<-gmat[,corr_ord]
        sink(paste(outpath,"/", prefix, ".", q_metadata.name,".Corr.xls",sep="")); cat("\t"); write.table(out,quote=FALSE,sep="\t");sink()
     }
   #------------------
   pdf(paste(outpath,"/", prefix, ".",colnames(q_metadata)[a],".Corr.pdf",sep=""),30,12)
   par(mfrow = c(4,10))
   for(m in 1:ncol(gmat)){
    par(mar = c(4,4,2,2))
    if(out[m, 2]<p.adj.cutoff && out[m, 1]>r.cutoff) {
      plot(q_metadata[,a], mat[,m], xlab=q_metadata.name[a], cex.lab=1.25, cex=4*out[m, 1],  ylab="Relative abundance",col=rgb(1,0,0,0.5), pch=20, main=colnames(mat)[m])
      mtext(paste("rho=",round(out[m, 1],2),"\n p.adj=",round(out[m, 2],2),sep=""),outer=FALSE,line=-3,adj=0.9,col="red")
      abline(glm(mat[,m]~q_metadata[,a]))
     }else if(out[m, 2]<p.adj.cutoff && out[m, 1]<(-r.cutoff)){
      plot(q_metadata[,a], mat[,m], xlab=q_metadata.name[a], cex.lab=1.25,cex=4*(-out[m, 1]),  ylab="Relative abundance",col=rgb(0,0.8,0.4,0.5), pch=20, main=colnames(mat)[m])
      mtext(paste("rho=",round(out[m, 1],2),"\n p.adj=",round(out[m, 2],2),sep=""),outer=FALSE,line=-3,adj=0.9,col="springgreen3")
      abline(glm(mat[,m]~q_metadata[,a]))
     }else{
      plot(q_metadata[,a], mat[,m], xlab=q_metadata.name[a], cex.lab=1.25, ylab="Relative abundance",col=rgb(0.6,0.6,0.6,0.5),pch=20, main=colnames(mat)[m])
      mtext(paste("rho=",round(out[m, 1],2),"\n p.adj=",round(out[m, 2],2),sep=""),outer=FALSE,line=-3,adj=0.9,col="grey60")
     }
     }
   dev.off()
}


}
