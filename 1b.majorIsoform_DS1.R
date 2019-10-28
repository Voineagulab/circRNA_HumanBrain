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

################ Major isoform analysis with random sampling
# generate 100 random samplings of circ and sj reads for each sample and calculate the circRNA and SJ expression as the average of the 100 random samplings
realC=data_DS1$circData_DS1$circCounts_filter[, -infocolsC]
realS=data_DS1$sjData_DS1$sjmaxCounts_filter[, -infocolsS]

# Note: "Sample" randomly samples from "Real"
sampleC=data_DS1$circData_DS1$circCounts_filter[, -infocolsC]
sampleS=data_DS1$sjData_DS1$sjmaxCounts_filter[, -infocolsS]
sampleC[, -infocolsC]=NA
sampleS[, -infocolsS]=NA
for (j in c(1:ncol(realC)))
{ 
sc=matrix(NA, ncol=1000, nrow=nrow(realC))
ss=matrix(NA, ncol=1000, nrow=nrow(realS))
for (i in c(1:1000))
{sc[,i]=rnegbin(n=nrow(realC),
                mu=mean(realC[,j]) ,
                theta=mean(realC[,j])^2/ (var (realC[,j]) - mean(realC[,j])))
ss[,i]=rnegbin(n=nrow(realS),
               mu=mean(realS[,j]) ,
               theta=mean(realS[,j])^2/ (var (realS[,j]) - mean(realS[,j])))
}
sampleC[,j]=apply(sc,1,mean)
sampleS[,j]=apply(ss,1,mean)
}
# calculate CLR for the random sampling data
# set Inf values to 1 and NA values to 0
clr_sample=sampleC; clr_sample[,]=NA
for (j in c(1:ncol(sampleC))) 
{
  clr_sample[,j]=sampleC[,j]/sampleS[,j]
  clr_sample[which(is.na(clr_sample[,j])==TRUE),j]=0
  clr_sample[which(clr_sample[,j]==Inf),j]=1
}
clr_sample=data.frame(clr_ds1[, infocolsC], clr_sample)

#############################Calculate major isoform relative clr
######Random sampling data
# calculate mean CLR across CTl samples for the random sampling and real data
clr_sample$MeanCTL=apply(clr_sample[, rownames(samples1[which(samples1$ASD.CTL=="CTL"),])] ,1,mean)
na1=grep("NA", clr_sample$EnsID)
clr_sample_agMax=aggregate(clr_sample$MeanCTL[-na1] ~ clr_sample$EnsID[-na1], FUN = max)
clr_sample_agSum=aggregate(clr_sample$MeanCTL[-na1] ~ clr_sample$EnsID[-na1], FUN = sum)
# Match the order of circRNAs between the CLR dataframes and the nclr_per_gene vector
nclr_per_gene_s=table(clr_sample$EnsID[-na1])
clr_sample_agMax=clr_sample_agMax[match(names(nclr_per_gene_s), clr_sample_agMax[,1]), ]
clr_sample_agSum=clr_sample_agSum[match(names(nclr_per_gene_s), clr_sample_agSum[,1]), ]
# Calculate major isoform relative CLR
major_isoform_ratio_clr_sample=clr_sample_agMax[,2]/clr_sample_agSum[,2]

########Real data
for (j in c((1+length(infocolsC)):ncol(clr_ds1))) 
{
  clr_ds1[which(is.na(clr_ds1[,j])==TRUE),j]=0
  clr_ds1[which(clr_ds1[,j]==Inf),j]=1
}
clr_ds1$MeanCTL=apply(clr_ds1[, rownames(samples1[which(samples1$ASD.CTL=="CTL"),])] ,1,mean)
na1=grep("NA", clr_ds1$EnsID)
clr_ds1_agMax=aggregate(clr_ds1$MeanCTL[-na1] ~ clr_ds1$EnsID[-na1], FUN = max)
clr_ds1_agSum=aggregate(clr_ds1$MeanCTL[-na1] ~ clr_ds1$EnsID[-na1], FUN = sum)
# Match the order of circRNAs in the CLR dataframes, and the nclr_per_gene vector
nclr_per_gene=table(clr_ds1$EnsID[-na1])
clr_ds1_agMax=clr_ds1_agMax[match(names(nclr_per_gene), clr_ds1_agMax[,1]), ]
clr_ds1_agSum=clr_ds1_agSum[match(names(nclr_per_gene), clr_ds1_agSum[,1]), ]
# Calculate major isoform relative CLR
major_isoform_ratio_clr_ds1=clr_ds1_agMax[,2]/clr_ds1_agSum[,2]

############Plot 
real=data.frame( nclr_per_gene,major_isoform_ratio_clr_ds1, "real" )
rownames(real)=real[,1]; real=real[,-1]
colnames(real)=c( "nIso", "ratio","categ")

sample=data.frame( nclr_per_gene_s,major_isoform_ratio_clr_sample, "sample" )
rownames(sample)=sample[,1]; sample=sample[,-1]
colnames(sample)=c( "nIso", "ratio","categ")

plotdata=as.data.frame(rbind(real, sample))
plotdata$categ=ifelse (plotdata$categ%in%"real", "DS1", "Random Sampling")

pdf("Analysis/Fig3/Fig3.MajorIsoform_with_RandomSampling_DS1_1000.pdf",height=4, width=10 )

ggplot(plotdata[which(plotdata$nIso < 26),], aes(x=factor(nIso), y=ratio, fill=categ)) +
  geom_boxplot() + 
  xlab("Number of circRNA isoforms") + 
  ylab("Major Isoform relative CLR") +
  theme(text=element_text( size=12))

dev.off()
############ Stats
for (i in sort(unique(nclr_per_gene)))
{
  use=which(nclr_per_gene==i)
  use_s=which(nclr_per_gene_s==i)
  print(i)
  print(wilcox.test(major_isoform_ratio_clr_ds1[use], major_isoform_ratio_clr_sample[use_s]))
  print("________")
}
