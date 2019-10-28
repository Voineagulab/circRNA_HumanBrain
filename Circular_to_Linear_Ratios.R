rm(list=ls())
library(WGCNA)
library(ggplot2);library(ggpubr)
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
###################################################################################################
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
par(pch=20, col="grey", family="serif")
############################################ LOAD DATA ############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/deconvolution_DS1.rda");load("DATA_TABLES/deconvolution_DS2.rda")
#Extract filtered data for genes, circRNAs and splice junctions,  DS1
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
circCts_ds1=data_DS1$circData_DS1$circCounts_filter; 
sjMCts_ds1=data_DS1$sjData_DS1$sjmaxCounts_filter
sjSCts_ds1=data_DS1$sjData_DS1$sjsumCounts_filter
#Extract filtered data for genes, circRNAs and splice junctions,  DS2
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
circCts_ds2=data_DS2$circData_DS2$circCounts_filter; 
sjMCts_ds2=data_DS2$sjData_DS2$sjmaxCounts_filter
sjSCts_ds2=data_DS2$sjData_DS2$sjsumCounts_filter
#infocols
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)
############################################ CALCULATE CIRCULAR TO LINEAR RATIOS #################
#Ratios are calculated as ratios of counts
# DS1
clr_ds1=circCts_ds1; ci_ds1=circCts_ds1
clr_ds1[,-infocolsC]=NA ; ci_ds1[,-infocolsC]=NA
for (j in c((length(infocolsC)+1): ncol(clr_ds1))) clr_ds1[,j]=circCts_ds1[,j]/sjMCts_ds1[,j]
for (j in c((length(infocolsC)+1): ncol(ci_ds1))) ci_ds1[,j]=circCts_ds1[,j]/(sjSCts_ds1[,j]/2+circCts_ds1[,j])
# DS2
clr_ds2=circCts_ds2; ci_ds2=circCts_ds2
clr_ds2[,-infocolsC]=NA ; ci_ds2[,-infocolsC]=NA
for (j in c((length(infocolsC)+1): ncol(clr_ds2))) clr_ds2[,j]=circCts_ds2[,j]/sjMCts_ds2[,j]
for (j in c((length(infocolsC)+1): ncol(ci_ds2))) ci_ds2[,j]=circCts_ds2[,j]/(sjSCts_ds2[,j]/2+circCts_ds2[,j])

#Save
ratios_ds1=list(clr=clr_ds1, ci=ci_ds1)
ratios_ds2=list(clr=clr_ds2, ci=ci_ds2)

save(ratios_ds1, file="DATA_TABLES/ratios_ds1.rda")
save(ratios_ds2, file="DATA_TABLES/ratios_ds2.rda")
