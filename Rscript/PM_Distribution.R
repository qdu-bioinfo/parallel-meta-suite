#!/usr/bin/env Rscript
#################################################################
# Function: Abundance distribution analysis
# Call: Rscript PM_Distribution.R -i abund_file -o outfile
# R packages used: optparse, reshape2, ggplot2,RColorBrewer, grDevices
# Authors: Yuzhu Chen, Zheng Sun, Xiaoquan Su, Honglei Wang 
# Updated at Aug. 20, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
# data_matrix -> Unclassified
#################################################################
# install necessary libraries
p <- c("optparse","reshape2","ggplot2","RColorBrewer")
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
  make_option(c("-i", "--abund_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),
  make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Optional]"),
  make_option(c("-o", "--out_dir"), type="character", default='Distribution', help="Output directory [default %default]"),
  make_option(c("-p", "--prefix"), type="character", default='Out', help="Output file prefix [default %default]"),
  make_option(c("-t", "--threshold"), type="double", default=0.01, help="Average value threshold [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$abund_file)) stop('Please input a feature table (*.Abd)')
# load data
matrixfile <- opts$abund_file
mapfile <- opts$meta_data
ave_t <- opts$threshold
outpath <- opts$out_dir
dir.create(outpath)
#------------------------------------------------------------------------------------
disbar <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
#-----------Edit
disbarm <- t(disbar)
#---------------
disbar <-disbar[order(rownames(disbar)),]

disbar <- t(disbar)

disbar <- disbar[names(sort(rowSums(disbar),decreasing = T)),]
disbar <- floor(disbar*1000000)/1000000
 Unclassified_other <- disbar[which(apply(disbar,1,mean) <= ave_t),]

invisible(if (sum( Unclassified_other) ==0 ) ( Unclassified_big <- disbar))
invisible(if (sum( Unclassified_other) !=0 ) ( Unclassified_big <- disbar[-(which(apply(disbar,1,mean) <= ave_t)),]))

widforpdf <- ncol(disbar)
 Unclassified_other <- as.matrix( Unclassified_other)
if (dim( Unclassified_other)[2] ==1 )  Unclassified_other <- t( Unclassified_other)

if (is.null( Unclassified_other)==F) {
  disbar <- rbind( Unclassified_big,"Other"=c(colSums( Unclassified_other)),deparse.level = 2)
}

if (mean(colSums(disbar))<0.9999) {                         #Complete to 100%
  if (rownames(disbar)[nrow(disbar)]=="Other") {
    disbar[nrow(disbar),] <- disbar[nrow(disbar),]+(1-colSums(disbar))
  }
  else {
    disbar <- rbind(disbar,"other"=(sapply((1-colSums(disbar)),function(x)max(x,0))),deparse.level = 2)
  }
}
#colours <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),sample(rainbow(length(colnames(t(disbar)))),length(colnames(t(disbar)))))
colours <- c(brewer.pal(4,"Set3")[1],brewer.pal(12,"Set3")[12],brewer.pal(12,"Set3")[3:6],brewer.pal(12,"Set3")[10:11],brewer.pal(8,"Set2")[2:8],brewer.pal(6,"Pastel2"),brewer.pal(8,"Pastel1"),brewer.pal(5,"Paired")[2:5],brewer.pal(3,"Accent"),sample(rainbow(length(colnames(t(disbar)))),length(colnames(t(disbar)))))

meltdata <- melt(abs(data.matrix(t(disbar))),varnames=c("Samples","Cutline"),value.name="Relative_Abundance")

#-----------------------------------------------------------------------------------------
if(is.null(opts$meta_data)) {
  pp<-ggplot(meltdata,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
    geom_bar(stat='identity')+ ylab("Relative Abundance")+ xlab("Samples")+
    scale_x_discrete(limits=c(colnames(disbarm)))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0%","25%","50%","75%","100%"))+
    guides(fill = guide_legend(ncol = (ceiling(nrow(disbar)/35))))+
    scale_fill_manual (values=colours) +
    theme(legend.position="right",axis.text.x=element_text(size=12,colour="black",angle=60,hjust=1,vjust=0.9),
          axis.text.y=element_text(size=12,colour="black"), axis.title.x=element_text(size=16,margin=margin(20,0,0,0)),
          axis.title.y=element_text(size=16,margin=margin(0,20,0,0)),panel.grid.major=element_line(colour=NA),legend.text = element_text(size = 10,colour = "black"),legend.title = element_text(size = 12,face="bold"),
		  legend.key.size = unit(0.7,'cm'),legend.key = element_rect(colour = "white", size=2),
		  plot.margin = unit(rep(1.5,4),'lines'),panel.background = element_blank(),axis.ticks.length = unit(0.25,'cm'))
  suppressWarnings(ggsave(paste(outpath, "/", opts$prefix, ".distribution.pdf",sep=""),plot=pp,width=ceiling(18+widforpdf/8),height=10, limitsize=FALSE))
}
if(is.null(opts$meta_data)==F) {
	data_map <- read.table(mapfile,header = T, row.names= 1,sep="\t",as.is=FALSE)
	data_map <- data_map[order(rownames(data_map)),]

	# Print OTU results
	pp<-ggplot(meltdata,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
		geom_bar(stat='identity')+ ylab("Relative Abundance")+ xlab("Samples")+
		scale_x_discrete(limits=c(colnames(disbarm)))+
		scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0%","25%","50%","75%","100%"))+
		guides(fill = guide_legend(ncol = (ceiling(nrow(disbar)/35))))+
		scale_fill_manual (values=colours) +
		theme(legend.position="right",
			axis.text.x=element_text(size=12,colour="black",angle=60,hjust=1,vjust=0.9),
			axis.text.y=element_text(size=12,colour="black"), 
			#margin: up right bottom left
			axis.title.x=element_text(size=16,margin=margin(20,0,0,0)),
			axis.title.y=element_text(size=16,margin=margin(0,20,0,0)),panel.grid.major=element_line(colour=NA),
			legend.text = element_text(size = 10,colour = "black"),legend.title = element_text(size = 12,face="bold"),
			plot.margin = unit(rep(1.5,4),'lines'),panel.background = element_blank(),
			legend.key.size = unit(0.7,'cm'),legend.key = element_rect(colour = "white", size=2),
			axis.ticks.length = unit(0.25,'cm'))
	suppressWarnings(ggsave(paste(outpath, "/", opts$prefix, ".distribution.pdf",sep=""),plot=pp,width=ceiling(18+widforpdf/8),height=10, limitsize=FALSE))
	########

	for (i in 1:ncol(data_map)) {
		if (is.factor(data_map[1,i])==T) {
			tempframe <- data.frame(meltdata,dv=data_map[,i])
			pp<-ggplot(tempframe,aes(x=Samples,y=Relative_Abundance,fill=Cutline))+
				geom_bar(stat = "identity")+ ylab("Relative Abundance")+ xlab("Samples")+
				#scale_x_discrete(limits=c(colnames(disbar)))+
				scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0%","25%","50%","75%","100%"))+
				guides(fill = guide_legend(ncol = (ceiling(nrow(disbar)/35))))+
				scale_fill_manual (values=colours) +
				facet_grid(~ dv,scales="free_x",space="free_x")+ 
				theme(legend.position="right",
					axis.text.x=element_text(size=12,colour="black",angle=60,hjust=1,vjust=0.95),
					axis.text.y=element_text(size=12,colour="black"), 
					axis.title.x=element_text(size=16,margin=margin(20,0,0,0)),
					axis.title.y=element_text(size=16,margin=margin(0,20,0,0)),panel.grid.major=element_line(colour=NA),
					strip.text.x = element_text(size=15),
					legend.key.size = unit(0.7,'cm'),legend.key = element_rect(colour = "white", size=2),
					legend.text = element_text(size = 10,colour = "black"),legend.title = element_text(size = 12,face="bold"),
					plot.margin = unit(rep(1.5,4),'lines'),axis.ticks.length = unit(0.25,'cm'))
			suppressWarnings(ggsave(paste(outpath,"/", opts$prefix, ".", colnames(data_map)[i],".distribution.pdf",sep=""),plot=pp,height=10,width=ceiling(18+widforpdf/8),limitsize=FALSE))
		}
	}
}


