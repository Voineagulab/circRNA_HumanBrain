rm(list=ls())
library(WGCNA)
library(ggplot2);library(ggpubr)
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION")
###################################################################################################
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
par(pch=20, col="grey", family="serif")
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)
############################################ LOAD DATA ############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/deconvolution_DS1.rda");load("DATA_TABLES/deconvolution_DS2.rda")
load("DATA_TABLES/ratios_ds1.rda");load("DATA_TABLES/ratios_ds2.rda")
de=read.csv("Analysis/DE/SuppTable.region_circ.csv")
r=de$Coordinate[grep("RIMS2", de$Symbol)]
ci1=t(ratios_ds1$ci[r,-infocolsC])
ci2=t(ratios_ds2$ci[r,-infocolsC])
pdf("Analysis/RIMS2_Boxplot.pdf", height=3.5, width=5.5)
colors=rep(c("lemonchiffon4","lemonchiffon1"), 10)
for(j in c(1:ncol(ci1)))
{
if (j==1) plotdata=list(CB1=ci1[grep("verm", rownames(ci1)), j],
              CTX1=ci1[-grep("verm", rownames(ci1)), j],
              CB2=ci2[grep("verm", rownames(ci2)), j],
              CTX2=ci2[-grep("verm", rownames(ci2)), j])
else plotdata=c(plotdata,list(CB1=ci1[grep("verm", rownames(ci1)), j],
                                       CTX1=ci1[-grep("verm", rownames(ci1)), j],
                                       CB2=ci2[grep("verm", rownames(ci2)), j],
                                       CTX2=ci2[-grep("verm", rownames(ci2)), j]) )
}

boxplot(plotdata, ylim=c(0,1), col=colors)
dev.off()