#!/usr/bin/env Rscript
#################################################################
# Function:  Random Forest analysis based on sequence count or sequence abd
# Call: Rscript PM_Marker_RFscore.R -i taxa.abd -m metadata -o outfile
# R packages used: randomForest ggplot2
# Authors: Yuzhu Chen, Zheng Sun, Xiaoquan Su
# Updated at Dec. 10, 2021
# Updated by Yuzhu Chen
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################
# install necessary libraries
p <- c("randomForest","ggplot2","optparse","RColorBrewer")
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
  make_option(c("-i", "--table_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),
  make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
  make_option(c("-o", "--out_dir"), type="character", default='RFimportance', help="Output file name[default %default]"),
  make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
  make_option(c("-n", "--ntree"), type="integer",default=5000, help="Ntree for Random Forest [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$table_file)) stop('Please input a feature (*.Abd)')
if(is.null(opts$meta_data)) stop('Please supply a meta data file')

# load data
matrixfile <- opts$table_file
mapfile <- opts$meta_data
outpath <- opts$out_dir

dir.create(outpath)

#------------------------------------------------------------------------------------
data_matrix <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
data_map <- read.table(mapfile,header = T, sep="\t", as.is=FALSE)

if (min(colSums(data_matrix))==0) {
  data_matrix <- data_matrix[,-which(apply(data_matrix,2,sum)==0)]              #trim 0 cols
}

setwd(outpath)
#------------------------------------------------------------------------------------
for (i in 2:ncol(data_map)) {
	if (is.numeric(data_map[1,i])==F) {
		tempframe=data.frame(data_matrix,status=as.factor(data_map[,i])) 
		#to regenerate a dataframe with discrete variable as last col                                                       
		bm.rf <- randomForest(status ~ .,data=tempframe,importance=T,proximity=T,ntree=opts$ntree)
		framebm <- data.frame(importance(bm.rf,type = 1))
		error_rate <- mean(bm.rf$confusion[,dim(bm.rf$confusion)[2]])
		IL <- ggplot(framebm,aes(y=reorder(row.names(framebm),framebm$MeanDecreaseAccuracy),x=framebm$MeanDecreaseAccuracy))+geom_point(size=3,colour="#D55E00")+geom_segment(aes(yend=row.names(framebm)),xend=0,colour="#D55E00")+theme_bw()+theme(panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())+xlab(paste("Importance (mean decrease in accuracy)")) + ylab(paste("Features"))+ggtitle(paste("Random-Forest Score (error rate = ",round(error_rate, digits = 4)*100,"%)",sep="")) + theme(title=element_text(size=10),plot.title=element_text(hjust = 0.5),axis.title.x=element_text(margin=margin(15,0,0,0)),axis.title.y=element_text(margin=margin(0,15,0,0)),plot.margin = unit(c(1,1.5,0.7,0.7),'lines'))
		newframebm <- data.frame(framebm,sn=row.names(framebm),sl=framebm[,1])
		outtxtfile <- paste(opts$prefix,".",names(data_map[i]),".RFimportance.txt", sep="")  #outputtxt
		cat(paste("Random-Forest Score (error rate = ",round(error_rate,digits = 4)*100,"%)\n",sep=""),file=outtxtfile)
		write.table(newframebm[order(newframebm[,1],decreasing=T),][,2:3],file=outtxtfile,quote=F,row.names = F,col.names = F,append = T)
		#outMDSfile <- paste(opts$prefix,"_",names(data_map[i]),".RFimportance_MDS.pdf", sep="")  #outputmds
		#pdf(outMDSfile)
		#suppressWarnings(MDSplot(bm.rf,tempframe$status))
		#invisible(dev.off())
		outPDFfile <- paste(opts$prefix,".",names(data_map[i]),".RFimportance.pdf", sep="")  #outputpdf
		suppressMessages(ggsave(filename=outPDFfile,plot=IL,width=9 ,height = (dim(newframebm)[1]*0.2 + 1),limitsize =F))
	} 
}
#invisible(file.remove("Rplots.pdf"))
invisible(system("rm -rf Rplots.pdf"));


