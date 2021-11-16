#!/usr/bin/env Rscript
#################################################################
# Function:  Multivariate statistical analysis based on sequence count or sequence abd
# Call: Rscript PM_ADiversity.R -i taxa.abd -m metadata -o outfile
# R packages used: optparse,vegan fossil abind
# Authors: Yuzhu Chen, Zheng Sun, Xiaoquan Su
# Updated at Aug. 20, 2021
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

# install necessary libraries
p <- c("optparse","vegan","fossil","abind","ggplot2","RColorBrewer")
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
	make_option(c("-i", "--table_file"), type="character", help="Input feature table with read count (*.Count) or abundance (*.Abd) [Required]"),
	make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
	make_option(c("-o", "--out_dir"), type="character", default='Alpha_diversity', help="Output directory [default %default]"),
	make_option(c("-w", "--width"), type="integer", default=10, help="Width of figure [default %default]"),
	make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$table_file)) stop('Please input a feature table (*.Count or *.Abd)')
if(is.null(opts$meta_data)) stop('Please input a meta data table')

# load data
matrixfile <- opts$table_file
mapfile <- opts$meta_data
outpath <- opts$out_dir

#------------------------------------------------------------------------------------
data_matrix <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
data_map <- read.table(mapfile,header = T, sep="\t",as.is=FALSE)

if (min(colSums(data_matrix))==0) {
	data_matrix <- data_matrix[,-which(apply(data_matrix,2,sum)==0)]              #trim 0 cols
}

ssn <- min(rowSums(data_matrix)) # standard sample number
indexes <- 0
chaos <- 0

#chao1
#for (i in 1:dim(data_matrix)[1])
#    chaos[i] <- chao1(data_matrix[i,])   # chao1 caculation
chaos <- apply(data_matrix, 1, chao1)

#shannon & simpson
indexes <- t(abind(diversity(data_matrix,index="shannon"),diversity(data_matrix,index="simpson"),chaos,along=0)) 

dir.create(outpath)
setwd(outpath)

outtxtfile <- paste(opts$prefix, ".Alpha_diversity_Index.txt", sep="")
outpdffile <- paste(opts$prefix, ".Alpha_diversity_Boxplot.pdf", sep="")

outputindex <- paste(rownames(data_matrix),"\t",indexes[,1],"\t",indexes[,2],"\t",indexes[,3])
write.table("Sample\tShannon\tSimpson\tChao1",file=outtxtfile,quote=F,,row.names=F,col.names = F)
write.table(outputindex,file=outtxtfile,append=T,quote=F,row.names=F,col.names = F)  #output

#-------------------------------boxplot-----------------------------------------
j=0
f= opts$width
x=c()
for (i in 2:ncol(data_map)) {
	if (is.null(data_map[1,i])==F) {
		f=nlevels(data_map[,i])+f
		if (nlevels(data_map[,i])!=0) (x <-c(x,nlevels(data_map[,i])))
		j=j+1
		if (nlevels(data_map[,i])==0) (x <-c(x,3))
	}
}  
colours <- c(brewer.pal(12,"Set3")[1:6],brewer.pal(12,"Set3")[10:11],brewer.pal(8,"Set2")[2:8],brewer.pal(6,"Pastel2"),brewer.pal(8,"Pastel1"),brewer.pal(5,"Paired")[2:5],brewer.pal(3,"Accent"),sample(rainbow(length(colnames(data_map))),length(colnames(data_map))))

dark_col <- c()
for (i in colours){
	rgbcol=col2rgb(i)
	red=rgbcol[1,1]
	green=rgbcol[2,1]
	blue=rgbcol[3,1]
	if(red-50<0){red = 0}else{red = red - 50}
	if(green-50<0){green = 0}else{green = green - 50}
	if(blue-50<0){blue = 0}else{blue = blue -50}
	dark=rgb(red,green,blue,maxColorValue=255)
	dark_col=c(dark_col,dark)
}


outpairwisefile <- paste(opts$prefix, ".Alpha_diversity_pairwised.txt", sep="")

write.table("Alpha diversity wilcoxon test\nPairwise.wilcox.test FDR",file=outpairwisefile,quote=F,,row.names=F,col.names = F)
pdf(outpdffile,height = 11,width=f*0.3*j)
layout(matrix(1:(j*3),3,,byrow = T),widths=x)
par(oma=c(1,2,1,2))
#shannon start
invisible(write.table("\n=====================\nShannon:",file=outpairwisefile,quote=F,,row.names=F,col.names = F,append = T,eol="\n\t"))
for (i in 2:ncol(data_map)) {
	if (is.factor(data_map[1,i])==T) {
		tempframe_shannon=data.frame(shan=c(indexes[,1]),dv=as.factor(data_map[,i])) 
		#to regenerate a dataframe with discrete variable as last col                                                       
		attach(tempframe_shannon)
		suppressWarnings(invisible(shant <- pairwise.wilcox.test(shan,dv,p.adjust="fdr")))
		suppressWarnings(write.table(shant$p.value,file=outpairwisefile,quote=F,row.names=T,col.names = T,append = T,sep = "\t"))
		suppressWarnings(write.table("",file=outpairwisefile,quote=F,row.names=F,col.names = F,append = T,eol = "\t"))
		invisible(if (nlevels(dv)==2) (mt=shant$p.value) else (mt=(kruskal.test(shan, dv))$p.value))
		boxplot(shan ~ dv, data = tempframe_shannon, xlab= colnames(data_map)[i],ylab="Shannon",main=paste("P=",round(mt,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2,col=colours,border=dark_col,medlwd=2,boxlwd=2,staplelwd=2)
		#                                    ggplot(tempframe_shannon,aes(y=shan,x=dv,fill=dv))+geom_boxplot()+theme_bw()+guides(fill=FALSE)+annotate("text",x=-Inf,y=Inf,label=paste("P=",round(mt,digits=4)),hjust=-.2,vjust=2)+xlab(colnames(data_map)[i])+ylab("Shannon")
		detach(tempframe_shannon)
	} 
	if (is.factor(data_map[1,i])==F){
		tempframe_shannon=data.frame(shan=c(indexes[,1]),dv=data_map[,i])
		attach(tempframe_shannon)
		suppressWarnings(corr <- cor.test(c(shan),c(dv),method = "pearson"))
		plot(shan,dv,pch=19,col=rgb(0,0,100,50,maxColorValue=255),ylab= colnames(data_map)[i],xlab="Shannon",main=paste("R=",round(corr$estimate,digits=4)," p=",round(corr$p.value,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2)
		lines(lowess(shan,dv),lwd=3,col="firebrick4")
		detach(tempframe_shannon)
	}
}
#simpson start
invisible(write.table("\n=====================\nSimpson:",file=outpairwisefile,quote=F,,row.names=F,col.names = F,append = T,eol="\n\t"))
for (i in 2:ncol(data_map)) {
	if (is.factor(data_map[1,i])==T) {
		tempframe_simpson=data.frame(shan=c(indexes[,2]),dv=as.factor(data_map[,i])) 
		#to regenerate a dataframe with discrete variable as last col              
		attach(tempframe_simpson)
		suppressWarnings(invisible(simpt <- pairwise.wilcox.test(shan,dv,p.adjust="fdr")))
		suppressWarnings(write.table(simpt$p.value,file=outpairwisefile,quote=F,row.names=T,col.names = T,append = T,sep = "\t"))
		suppressWarnings(write.table("",file=outpairwisefile,quote=F,row.names=F,col.names = F,append = T,eol = "\t"))
		invisible(if (nlevels(dv)==2) (mt=simpt$p.value) else (mt=(kruskal.test(shan, dv))$p.value))
		boxplot(shan ~ dv, data = tempframe_simpson, xlab= colnames(data_map)[i],ylab="Simpson",main=paste("P=",round(mt,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2,col=colours,border=dark_col,medlwd=2,boxlwd=2,staplelwd=2)
		#                                    ggplot(tempframe_simpson,aes(y=simpt,x=dv,fill=dv))+geom_boxplot()+theme_bw()+guides(fill=FALSE)+annotate("text",x=-Inf,y=Inf,label=paste("P=",round(mt,digits=4)),hjust=-.2,vjust=2)+xlab(colnames(data_map)[2])+ylab("Simpson")
		detach(tempframe_simpson)
	}
	if (is.factor(data_map[1,i])==F){
		tempframe_simpson=data.frame(shan=c(indexes[,2]),dv=data_map[,i])
		attach(tempframe_simpson)
		suppressWarnings(corr <- cor.test(c(shan),c(dv),method = "pearson"))
		plot(shan,dv,pch=19,col=rgb(0,0,100,50,maxColorValue=255),ylab= colnames(data_map)[i],xlab="Simpson",main=paste("R=",round(corr$estimate,digits=4)," p=",round(corr$p.value,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2)
		lines(lowess(shan,dv),lwd=3,col="deepskyblue4")
		detach(tempframe_simpson)
	}
}
#chao1 start
invisible(write.table("\n=====================\nCHAO1:",file=outpairwisefile,quote=F,,row.names=F,col.names = F,append = T,eol="\n\t"))
for (i in 2:ncol(data_map)) {
	if (is.factor(data_map[1,i])==T) {
		tempframe_chaoooo=data.frame(shan=c(indexes[,3]),dv=as.factor(data_map[,i])) #to regenerate a dataframe with discrete variable as last col        
		attach(tempframe_chaoooo)
		suppressWarnings(invisible(chaot <- pairwise.wilcox.test(shan,dv,p.adjust="fdr")))
		suppressWarnings(write.table(chaot$p.value,file=outpairwisefile,quote=F,row.names=T,col.names = T,append = T,sep = "\t"))
		suppressWarnings(write.table("",file=outpairwisefile,quote=F,row.names=F,col.names = F,append = T,eol = "\t"))
		invisible(if (nlevels(dv)==2) (mt=chaot$p.value) else (mt=(kruskal.test(shan, dv))$p.value))
		boxplot(shan ~ dv, data = tempframe_chaoooo, xlab= colnames(data_map)[i],ylab="CHAO1",main=paste("P=",round(mt,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2,col=colours,border=dark_col,medlwd=2,boxlwd=2,staplelwd=2)
		#                                    ggplot(tempframe_chaoooo,aes(y=chaot,x=dv,fill=dv))+geom_boxplot()+theme_bw()+guides(fill=FALSE)+annotate("text",x=-Inf,y=Inf,label=paste("P=",round(mt,digits=4)),hjust=-.2,vjust=2)+xlab(colnames(data_map)[2])+ylab("CHAO1")
		detach(tempframe_chaoooo)
	}
	if (is.factor(data_map[1,i])==F){
		tempframe_chaoooo=data.frame(shan=c(indexes[,3]),dv=data_map[,i])
		attach(tempframe_chaoooo)
		suppressWarnings(corr <- cor.test(c(shan),c(dv),method = "spearman"))
		plot(shan,dv,pch=19,col=rgb(0,0,100,50,maxColorValue=255),ylab= colnames(data_map)[i],xlab="Chao",main=paste("R=",round(corr$estimate,digits=4)," p=",round(corr$p.value,digits=4)),cex.main=1.4,cex.lab=1.3,cex.axis=1.2)
		lines(lowess(shan,dv),lwd=3,col="red3")
		detach(tempframe_chaoooo)
	}
}
invisible(dev.off())
