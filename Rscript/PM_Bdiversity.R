#!/usr/bin/env Rscript
#################################################################
# Function: Multivariate statistical analysis based on distance matrix
# Call: Rscript PM_Bdiversity.R -m map_file -d dist_file -o output
# R packages used: reshape,ggplot2,pheatmap,pROC,combinat,plyr,vegan,optparse,parallel
# Authors: Yuzhu Chen, Mingqian Zhang
# Updated at Aug. 20, 2021
# Updated by Yuzhu Chen, Mingqian Zhang
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

options(warn=-1)
#Rprof()
## install necessary libraries
p <- c("reshape","ggplot2","pheatmap","pROC","combinat","plyr","vegan","optparse","parallel")

usePackage <- function(p){
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
sourcedir <- Sys.getenv("ParallelMETA")

c_thread<-detectCores()

# make option list and parse command line
option_list <- list(
	make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix table [Required]"),
	make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
	make_option(c("-o", "--out_dir"), type="character", default='Beta_diversity', help="Output directory [default %default]"),
	make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
	make_option(c("-n", "--dist_name"), type="character", default='Default', help="The distance metrics name such as Meta-Storms, Jensen-Shannon, Euclidean et al. [Optional, default %default]"),    
	make_option(c("-t", "--threads"),type="numeric",default='c_thread',help="Number of thread [default is auto]")#zmq
)

opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$meta_data)) stop('Please input a meta data file')
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

# create output directory
dir.create(outpath1<-paste(opts$out_dir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

# read arguments
dist_matrix_name<-opts$dist_file                       
meta_name<-opts$meta_data               
dist_name<-opts$dist_name
prefix<-opts$prefix 

if(is.na(opts$threads)){
	opts$threads<-c_thread
}

t<-opts$threads

# read tables
dist_matrix <- read.table(dist_matrix_name,header=T,row.names=1)
meta_data <- read.table(meta_name,header=T,sep="\t",row.names=1,as.is=FALSE)

#write anosim adonis process into xxx.Beta_diversity_log.txt
con <- file(paste(opts$out_dir,'/',prefix,'.Beta_diversity_log.txt',sep=''))
sink(con, append=TRUE)
sink(con, append=TRUE, type='message')

if(length(meta_data) == 1){
	all_group <- colnames(meta_data)
}else{
	all_group <- colnames(meta_data)
	all_group_f<-colnames(meta_data)[sapply(meta_data,class)=="factor"]
	all_group_n<-colnames(meta_data)[sapply(meta_data,class)!="factor"]
}
cat("All the sample metadata: ",all_group, "\n\n",sep=" ")

# Data Check
if(any((colnames(dist_matrix)==rownames(dist_matrix))==FALSE)) 
	cat("The column names do not exactly match the row names! Please revise!")
if(any((rownames(meta_data)==rownames(dist_matrix))==FALSE)) 
	cat("The row names in Map file do not exactly match the row names in the distance matrix! Please revise!\n")
	

###### Anosim and Adonis test
stat_summ<-matrix(NA,nrow=length(all_group),ncol=4)
rownames(stat_summ)<-all_group;
colnames(stat_summ)<-c("Adonis.F","Adonis.P","Anosim.R","Anosim.P")
#--------------------------------
#suppressWarnings(
for(group in all_group) {
	#--------------------------------
	if(is.element(group,all_group_f)){
		ano<-anosim(dist_matrix,meta_data[,group],parallel=t)
		stat_summ[group,4]<-ano.P<-ano$signif
		stat_summ[group,3]<-ano.R<-ano$statistic
		cat("ANOSIM (",group,"): \n")
		cat("--------------------------------")
		print(ano)
	}
	#--------------------------------
	ado<-adonis(dist_matrix~meta_data[,group],parallel=t)
	stat_summ[group,2]<-ado.P<-ado$aov.tab$P[1]
	stat_summ[group,1]<-ado.F<-ado$aov.tab$F.Model[1]
	cat("ADONIS/PERMANOVA (",group,"): \n")
	cat("--------------------------------\n")
	print(ado$aov.tab)
	cat("--------------------------------\n\n")
}
#)
#write anosim adonis results into xxxx.Beta_diversity_summ.xls (a table include the result)
sink(paste(outpath1,prefix,".Beta_diversity_summ.xls",sep=""));cat("\t");write.table(stat_summ,quote=FALSE,sep='\t',row.names=TRUE);sink()
###### two test finish

##

#--------------------------------
# Categorical metadata: Distance boxplot
#--------------------------------
if(length(all_group_f)>=1){
	for(g in all_group_f){
		dir.create(outpath2<-paste(outpath1,prefix,".",g,"/",sep=""))
		dm<-dist_matrix
		group<-meta_data[,g]
		group_name=g
		## data check
		if(ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE)
			stop('The distance matrix is not squared')
		if( length(unique(group))==1)
			stop('At least two levels for a given sample category in your metadata file are required.')
		if( length(group)!=nrow(dm))
			stop('The number of rows in metadata and distance matrix are not equal')
		if(nlevels(group)>length(group)*0.9)
			stop('The number of levels in a certain category can not exceed 90% of total number of samples')
		dm<-dm[order(group),order(group)]
		mt<-group[order(group)]
		names(mt)<-rownames(dm)
		colnames(dm)<-rownames(dm)<-paste(rownames(dm),mt,sep="____") 
		dm[lower.tri(dm)]<-NA
		melt_dm<-melt(data.matrix(dm))
		rm(dm)
		melt_dm<-melt_dm[!is.na(melt_dm$value),]
		melt_dm<-melt_dm[which(melt_dm$X1!=melt_dm$X2),]
		Row_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm$X1),"____",fixed=TRUE)))
		Col_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm$X2),"____",fixed=TRUE)))
		VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
		DistType<-as.factor(Row_Info[,2]==Col_Info[,2])
		DistType<-factor(DistType,levels=levels(DistType),labels=c("AllBetween","AllWithin"))
		plot_value<-data.frame(DistType,VS,Row_Info,Col_Info,d=melt_dm$value)
		colnames(plot_value)<-c("DistType","GroupPair","Sample_1","Group_1","Sample_2","Group_2","Dist")
		#plot_value<-data.frame(VS,DistType,d=melt_dm$value)
		#colnames(plot_value)<-c("GroupPair","DistType","Dist")
		rm(Row_Info)
		rm(Col_Info)
		rm(VS)
		rm(DistType)
		rm(melt_dm)
		# Output distance table and boxplot
		## grouppair>4,two test 
		if(length(plot_value$Dist)>4){
			#t-test  p-value
			p_t<-with(plot_value,t.test(Dist~DistType))$p.value
			#Wilcoxon rank sum test  p-value
			p_w<-with(plot_value,wilcox.test(Dist~DistType))$p.value
		}else{
			p_t<-NA
			p_w<-NA
		}
		## Output distance-values table
		filepath<-sprintf('%s%s%s%s%s',outpath2,dist_name,'.',group_name,'.plot_values.xls')
		sink(filepath);write.table(plot_value,quote=FALSE,sep='\t',row.names=FALSE);sink()
		## Output distance boxplot
		if(nlevels(group)<30){
			plot<-qplot(x=GroupPair, y=Dist, data=plot_value, geom='boxplot',position='dodge',main='',xlab="Group pair",ylab=paste(dist_name,' Distance',sep=''),outlier.alpha=0) + coord_flip() +  theme_bw() + theme(axis.title.x=element_text(margin=margin(15,0,0,0)),axis.title.y=element_text(margin=margin(0,15,0,0)),plot.margin = unit(c(0.5,1.3,0.8,0.9),'lines'))
			suppressMessages(ggsave(filename=paste(outpath2,dist_name,'.',group_name,'.boxplot.ggplot.pdf',sep=''),plot=plot,height=ifelse(nlevels(mt)>2,nlevels(mt),2))) 
		}
		rm(plot)
		b_plot_value<-data.frame(Grouping=rep(g,nrow(plot_value)),plot_value)
		if(g==all_group_f[1]){
			bp_value<-b_plot_value
		}else{
			bp_value<-rbind(bp_value,b_plot_value)
		}
		rm(plot_value)
		rm(b_plot_value)
	}
	p<-qplot(x=Grouping, y=Dist, data=bp_value, geom="boxplot", fill=DistType, position="dodge",main="", ylab=paste(dist_name," Distance",sep=""),outlier.alpha=0)+coord_flip()+ theme_bw() +theme(axis.title.x=element_text(margin=margin(15,0,0,0)),axis.title.y=element_text(margin=margin(0,15,0,0)),panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),plot.margin = unit(rep(0.6,4),'lines'))
	suppressMessages(ggsave(filename=paste(outpath1,"/",prefix,".Beta_diversity_DistBoxplot.pdf",sep=""),plot=p, limitsize=TRUE, width=6, height=ifelse(length(all_group_f)>1,length(all_group_f)*1.4,2)))
	rm(p)
	rm(bp_value)
}


#--------------------------------
# Continuous metadata: Scatterplot
#--------------------------------
if(length(all_group_n)>=1){
	for(group in all_group_n){
		dir.create(outpath2<-paste(outpath1,prefix,".",group,"/",sep=""))
		
		dist_num<-data.matrix(dist(meta_data[,group]))
		dm_n_value<-dist_num[lower.tri(dist_num, diag = FALSE)]
		dm_value<-dist_matrix[lower.tri(dist_matrix, diag = FALSE)]
		
		dm_data<-data.frame(dm=dm_value,group=dm_n_value)
		
		corr<-with(dm_data,cor.test(dm_value,dm_n_value,method="spearman"))
		
		rm(dist_num)
		rm(dm_n_value)
		rm(dm_value)
		
		p<-ggplot(data=dm_data, aes(x=group,y=dm))+geom_point(alpha=0.2)+annotate("text", x=max(dm_data$group)*0.9, y=max(dm_data$dm)*0.9, label= paste("Rho=",round(corr$estimate,2),"\n","P=",round(corr$p.value,4),"\n",sep=""))+ ylab(paste(dist_name," Distance",sep="")) + xlab(bquote(paste(~Delta, .(group),sep=" ")))+theme_bw()+ theme(axis.title.x=element_text(margin=margin(10,0,0,0)),axis.title.y=element_text(margin=margin(0,10,0,0)),plot.margin = unit(c(1,1,0.5,0.5),"lines")) 
		
		if(corr$estimate>0.4 && corr$p.value<0.01) p<-p+geom_smooth(method = "loess", se=TRUE, span=1) 
		suppressMessages(ggsave(filename=paste(outpath2,"/",prefix,".Beta_diversity_",group,".Scatterplot",".pdf",sep=""),plot=p, limitsize=TRUE, width=4, height=4))
		
		rm(p)
	}
}

sink()
sink(type='message')

