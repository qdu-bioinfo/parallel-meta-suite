#!/usr/bin/env Rscript
#################################################################
# Function: Rarefaction curve
# Call: Rscript Rarefaction.R -i count_table -o outfile
# R packages used: optparse
# Last update: 2015-09-18, Zheng Sun, Xiaoquan Su
#################################################################

## install necessary libraries
p <- c("optparse","vegan","plyr","permute","lattice","reshape2","ggplot2","squash","RColorBrewer","scales","grid" )
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default='Rarefaction', help="Output path [default %default]"),
  make_option(c("-p", "--prefix"), type="character", default='Out',help="Output file prefix [Optional, default %default]"),
  make_option(c("-a", "--annotation"), type="logical", default='T',help="annotation [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

annotation <- opts$annotation
outpath <- opts$outpath
outpath <- paste(outpath,"/",sep="")
outprefix <- opts$prefix

## clean R environment
setwd(outpath)


data_matrix <- read.table(file="Shannon.txt",header = T,sep="\t")
data_max <- read.table(file="Shannon_max.txt",header = T,sep="\t")
numberofsamples <- dim(data_max)[1]
colors<-sample(rainbow(numberofsamples),numberofsamples)
po1 <- ggplot(data_matrix,aes(x=x,y=y))+geom_line(size=0.8,aes(colour=SampleID)) + theme_bw() +
  scale_colour_manual(values=colors) + 
  guides(colour=guide_legend(title="Samples",ncol=ceiling(numberofsamples/25))) +
  xlab("Number of Sequence") + ylab("Shannon Index") + ggtitle("Rarefaction curves") +
  theme(legend.text=element_text(size=10),legend.key.size=unit(1.2,"cm"),
        panel.border=element_rect(fill=NA),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=14,colour="black"),
        axis.text.y=element_text(size=14,colour="black"),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24)) 
if(annotation) {
po1 <- po1 + annotate("text",x=as.numeric(data_max[,2]),y=as.numeric(data_max[,3]),label=data_max[,1],size=3,hjust=-0.5)}
suppressWarnings(ggsave(filename =paste(outprefix,".Rarefaction_shannon.pdf",sep=""),plot = po1,width=(18+ceiling(numberofsamples/25)*2),height=15))

## for otu curve
data_matrix <- read.table(file="Observe_otu.txt",header = T,sep="\t")
data_max <- read.table(file="Observe_otu_max.txt",header = T,sep="\t")
numberofsamples <- dim(data_max)[1]
colors<-sample(rainbow(numberofsamples),numberofsamples)
po2 <- ggplot(data_matrix,aes(x=x,y=y))+geom_line(size=0.8,aes(colour=SampleID)) + theme_bw() +
  scale_colour_manual(values=colors) + 
  guides(colour=guide_legend(title="Samples",ncol=ceiling(numberofsamples/25))) +
  xlab("Number of Sequence") + ylab("Number of Observed Taxa") + ggtitle("Rarefaction curves") +
  theme(legend.text=element_text(size=10),legend.key.size=unit(1.2,"cm"),
        panel.border=element_rect(fill=NA),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=14,colour="black"),
        axis.text.y=element_text(size=14,colour="black"),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24)) 
if(annotation) {
po2 <-po2 + annotate("text",x=as.numeric(data_max[,2]),y=as.numeric(data_max[,3]),label=data_max[,1],size=3,hjust=-0.5)}
suppressWarnings(ggsave(filename =paste(outprefix,".Rarefaction.pdf",sep=""),plot = po2,width=(18+ceiling(numberofsamples/25)*2),height=15))
