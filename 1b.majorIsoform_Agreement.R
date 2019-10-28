rm(list=ls())
library(MASS)
library(data.table)
library(WGCNA)
library(ggplot2);
library(ggpubr); 
library(dplyr);

setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
###################################################################################################
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
par(pch=20, col="grey", family="serif")
############################################ LOAD DATA ############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/deconvolution_DS1.rda");load("DATA_TABLES/deconvolution_DS2.rda")
load("DATA_TABLES/ratios_ds1.rda");load("DATA_TABLES/ratios_ds2.rda")
#Extract filtered data for genes, circRNAs and splice junctions,  DS1
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
circCts_ds1=data_DS1$circData_DS1$circCounts_filter; circCpm_ds1=data_DS1$circData_DS1$circCpm_filter
geneCts_ds1=data_DS1$geneData_DS1$counts_filter; geneRPKM_ds1=data_DS1$geneData_DS1$rpkm_filter
sjMCpm_ds1=data_DS1$sjData_DS1$sjmaxCpm_filter;sjMCts_ds1=data_DS1$sjData_DS1$sjmaxCounts_filter
sjSCpm_ds1=data_DS1$sjData_DS1$sjsumCpm_filter;sjSCts_ds1=data_DS1$sjData_DS1$sjsumCounts_filter

samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
circCts_ds2=data_DS2$circData_DS2$circCounts_filter; circCpm_ds2=data_DS2$circData_DS2$circCpm_filter
geneCts_ds2=data_DS2$geneData_DS2$counts_filter; geneRPKM_ds2=data_DS2$geneData_DS2$rpkm_filter
sjMCpm_ds2=data_DS2$sjData_DS2$sjmaxCpm_filter;sjMCts_ds2=data_DS2$sjData_DS2$sjmaxCounts_filter
sjSCpm_ds2=data_DS2$sjData_DS2$sjsumCpm_filter;sjSCts_ds2=data_DS2$sjData_DS2$sjsumCounts_filter

clr_ds1=ratios_ds1$clr; ci_ds1=ratios_ds1$ci
clr_ds2=ratios_ds2$clr; ci_ds2=ratios_ds2$ci

# Define infocols
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)

##################Agreement
use_clr_ds1=clr_ds1[-grep("NA", clr_ds1$EnsID) ,]
n_per_gene1=as.data.frame(table(use_clr_ds1$EnsID))
isoform=n_per_gene1$Var1[(which(n_per_gene1$Freq <= 300) & (n_per_gene1$Freq > 1))]
use_clr_ds1=use_clr_ds1[which(use_clr_ds1$EnsID%in%isoform) ,]
#use_clr_ds1=use_clr_ds1[-which(use_clr_ds1$EnsID%in%isoform) ,]
#use_clr_ds1=use_clr_ds1[which(rowSums(use_clr_ds1[,-infocolsC] > 0) > 144/2 ),]

rank1=use_clr_ds1 ; rank1[,-infocolsC]=NA
for (j in c( (length(infocolsC)+1) : ncol(rank1)))
{
temp=use_clr_ds1[, c(infocolsC, j)]
k=ncol(temp)
temp[which(is.na(temp[,k]) ==TRUE),k]=0
temp[which(temp[,k] ==Inf ),k ]=1
clrMax=aggregate(temp[,k] ~ temp$EnsID, FUN = max); colnames(clrMax)=c("EnsID", "Max")
m=match(temp$EnsID, clrMax$EnsID)
rank1[,j]=ifelse(temp[,k] == clrMax$Max[m], 1, -1)
rank1[which(clrMax$Max[m]==0),j]=0
rm(temp)
}


use_clr_ds2=clr_ds2[-grep("NA", clr_ds2$EnsID) , ]
n_per_gene2=as.data.frame(table(use_clr_ds2$EnsID))
isoform=n_per_gene2$Var1[which((n_per_gene2$Freq <= 300))]
use_clr_ds2=use_clr_ds2[which(use_clr_ds2$EnsID%in%isoform) ,]
#use_clr_ds2=use_clr_ds2[-which(use_clr_ds2$EnsID%in%isoform) ,]
#use_clr_ds1=use_clr_ds1[which(rowSums(use_clr_ds1[,-infocolsC] > 0) > 144/2 ),]
rank2=use_clr_ds2 ; rank2[,-infocolsC]=NA
for (j in c( (length(infocolsC)+1) : ncol(rank2)))
{
  temp=use_clr_ds2[, c(infocolsC, j)]
  k=ncol(temp)
  temp[which(is.na(temp[,k]) ==TRUE),k]=0
  temp[which(temp[,k] ==Inf ),k ]=1
  clrMax=aggregate(temp[,k] ~ temp$EnsID, FUN = max); colnames(clrMax)=c("EnsID", "Max")
  m=match(temp$EnsID, clrMax$EnsID)
  rank2[,j]=ifelse(temp[,k] == clrMax$Max[m], 1, -1)
  rank2[which(clrMax$Max[m]==0),j]=0
  rm(temp)
}

mi1=which(rowSums(rank1==1)/(rowSums(rank1==1)+ rowSums(rank1==-1)) > 0.5)
m=match(rownames(rank1[mi1,]), rownames(rank2))

mi2=rowSums(rank2==1)/(rowSums(rank2==1) +rowSums(rank2==-1))
pdf("Analysis/SuppFigs/MajorIsoformCallAgreement.pdf", height=4,width=4)
hist(mi2[m], xlab="Proportion of DS2 samples", col="indianred", main="Major isoform agreement in DS1 vs. DS2")
abline(v=0.5, lty=2)
dev.off()
# table(mi2[m] > 0.5)
# table(is.na(m))
# hist(mi1)
# propMajor1=rep(0, nrow(rank1))
# for(j in c(1:nrow(rank1)))
#   {x=rank1[j,]
#   propMajor1[j]=length(which(x==1))/(length(x)-length(which(x==0)))}
# hist(propMajor1)
# boxplot(propMajor1  ~ n_per_gene1$Freq[match(use_clr_ds1$EnsID, n_per_gene1$Var1)])
# abline
# propMajor2=rep(0, nrow(rank2))
# for(j in c(1:nrow(rank2)))
# {x=rank2[j,]
# propMajor2[j]=length(which(x==1))/(length(x)-length(which(x==0)))}
# hist(propMajor2)
# #test=apply(t(rank1[, -infocolsC]), 2, as.data.frame(table))
# 
# perc_agreement=function(x, mx)
# {
# result=rep(NA, ncol(mx))
# for (j in c(1:length(result)))
#   {major=which(x==1)
#   result[j]=length(which(mx[major,j]==1)) / (length(which(mx[major,j]==1)) +length(which(mx[major,j]== -1)) )}
# return(result)
# }
# 
# use_rank1=rank1[,-infocolsC]
# ds1_perc_agreement=matrix(NA, ncol=ncol(use_rank1), nrow=ncol(use_rank1))
# colnames(ds1_perc_agreement)=rownames(ds1_perc_agreement)=colnames(use_rank1)
# for(j in c(1:ncol(use_rank1)))
#   ds1_perc_agreement[,j]=perc_agreement(use_rank1[,j], use_rank1)
# 
# use_rank2=rank2[,-infocolsC]
# ds2_perc_agreement=matrix(NA, ncol=ncol(use_rank2), nrow=ncol(use_rank2))
# colnames(ds2_perc_agreement)=rownames(ds2_perc_agreement)=colnames(use_rank2)
# for(j in c(1:ncol(use_rank2)))
#   ds2_perc_agreement[,j]=perc_agreement(use_rank2[,j], use_rank2)
# 
# boxplot(ds1_perc_agreement, ylim=c(0,1), col="red")
# boxplot(ds2_perc_agreement, ylim=c(0,1), col="blue")
# f=function(x)
# {
# t=table(t(x))
# main_iso=names(t)[which(t==max(t))]
# if (is.na(t["1"])==TRUE) perc_major =0 else
# {
# if (is.na(t["-1"])==FALSE) perc_major=t["1"]/(t["1"]+t["-1"]) else perc_major=1
# }
# return(perc_major)
# }
# total_rank1=apply(rank1[,-c(infocolsC, ncol(rank1))], 1, f)
# total_rank1=unlist(total_rank1)
# # total_rank1=transpose(total_rank1)[[1]]
# # total_rank1=as.numeric(total_rank1)
# # names(total_rank1)=rownames(rank1)
# 
# total_rank2=apply(rank2[,-infocolsC], 1, f)
# total_rank2=unlist(total_rank2)
# # total_rank2=transpose(total_rank2)[[1]]
# # total_rank2=as.numeric(total_rank2)
# # names(total_rank2)=rownames(rank2)
# 
# common=intersect(names(total_rank1) , names(total_rank2))
# total_rank1=total_rank1[match(common, names(total_rank1))]
# total_rank2=total_rank2[match(common, names(total_rank2))]
# 
# plotdata2=data.frame(total_rank1, total_rank2)
# ggplot(plotdata2, aes(x=total_rank1, y=total_rank2)) + geom_point(col="midnightblue") +geom_smooth()
