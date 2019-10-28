setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
rm(list=ls())
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
library(splines)
par(pch=20, col="grey", family="serif")
############################################ LOAD DATA
############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/ratios_ds1.rda");load("DATA_TABLES/ratios_ds2.rda")
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
infocolsG=c(1:2); infocolsC=c(1:15)

CON=as.character(samples1$Sample[which(samples1$ASD.CTL == "CTL")])
ASD=as.character(samples1$Sample[which(samples1$ASD.CTL == "ASD")])
# genes=data_DS1$geneData_DS1$rpkm_filter
# circ=data_DS1$circData_DS1$circCpm[rowSums(data_DS1$circData_DS1$circCpm[,CON] > 0) > 34 , ]
dataG=data_DS1$geneData_DS1$rpkm_filter; dataS=data_DS1$sjData_DS1$sjmaxCpm_filter; dataC=data_DS1$circData_DS1$circCpm_filter

meanG=apply(dataG[, CON], 1,mean);sdG=apply(dataG[, CON], 1,sd)
meanC=apply(dataC[, CON], 1,mean);sdC=apply(dataC[, CON], 1,sd)
meanS=apply(dataS[, CON], 1,mean);sdS=apply(dataS[, CON], 1,sd)

nSamples_G=rowSums(dataG[,CON] > 0)
nSamples_C=rowSums(dataC[,CON] > 0)
nSamples_S=rowSums(dataS[,CON] > 0)

# write.csv(meanG, "Analysis/Fig5/Fig5A.MeanG.csv")
# write.csv(meanS, "Analysis/Fig5/Fig5A.MeanS.csv")
# write.csv(meanC, "Analysis/Fig5/Fig5A.MeanC.csv")
# 
# write.csv(sdG, "Analysis/Fig5/Fig5A.sdG.csv")
# write.csv(sdS, "Analysis/Fig5/Fig5A.sdS.csv")
# write.csv(sdC, "Analysis/Fig5/Fig5A.sdC.csv")
# 
# write.csv(nSamples_G, "Analysis/Fig5/Fig5B.nSamples_G.csv")
# write.csv(nSamples_C, "Analysis/Fig5/Fig5B.nSamples_C.csv")
# write.csv(nSamples_S, "Analysis/Fig5/Fig5B.nSamples_S.csv")


pdf("Analysis/SuppFigs/SuppFig.Mean_vs_Sd.pdf", height=8, width=8)
par(pch=20, col="slategrey", family="serif")
par(mfrow=c(3,3))
cvC=sdC/meanC
cvG=sdG/meanG
cvS=sdS/meanS
plot(log2(sdG) ~ log2(meanG), sub=paste("CV=", round(mean(cvG, na.rm=TRUE),2))); abline(0,1, col="red")
plot(log2(sdS) ~ log2(meanS),sub=paste("CV=", round(mean(cvS, na.rm=TRUE),2))); abline(0,1, col="red" )
plot(log2(sdC) ~ log2(meanC),sub=paste("CV=", round(mean(cvC, na.rm=TRUE),2)), xlim=c(-12.5,5), ylim=c(-12.,5)); abline(0,1, col="red" )
wilcox.test(cvC, cvS)
dev.off()

pdf("Analysis/SuppFigs/SuppFig.nSamples_hist.pdf", height=8, width=8)
par(mfrow=c(3,3))
hist(nSamples_G, col="lightgrey")
hist(nSamples_S, col="lightgrey")
hist(nSamples_C, col="lightgrey")
dev.off()

pdf("Analysis/Fig5/Fig5.ScatterPlots.pdf", height=4, width=8)
par(mfrow=c(1,2))
circCpm_summ_ds1=calcSummaryData(data_DS1$circData_DS1$circCpm_filter, infocolsC, samples1, th=0.1)
circCpm_summ_ds2=calcSummaryData(data_DS2$circData_DS2$circCpm_filter, infocolsC, samples2, th=0.1)
clr_summ_ds1=calcSummaryData(ratios_ds1$clr, infocolsC, samples1, 0.01)
clr_summ_ds2=calcSummaryData(ratios_ds2$clr, infocolsC, samples2, 0.01)

m=match(rownames(circCpm_summ_ds1), rownames(circCpm_summ_ds2))
use=which(is.na(m)==FALSE)
cor(circCpm_summ_ds1$nexpPerFeature/nrow(samples1), circCpm_summ_ds2$nexpPerFeature[m]/nrow(samples2), use="pair", method="s")
plot(circCpm_summ_ds1$nexpPerFeature/nrow(samples1), circCpm_summ_ds2$nexpPerFeature[m]/nrow(samples2), pch=20,col="slategray3", xlab="Proportion of DS1 Samples", ylab="Proportion of DS2 Samples",  main=" ")
write.csv(cbind(circCpm_summ_ds1$Coordinate[use], circCpm_summ_ds1$nexpPerFeature[use]/nrow(samples1), circCpm_summ_ds2$nexpPerFeature[m[use]]/nrow(samples2)) , "SourceData/Fig5C.csv")
m=match(rownames(clr_summ_ds1), rownames(clr_summ_ds2))
cor(clr_summ_ds1$Mean, clr_summ_ds2$Mean[m], use="pair", method="s")
plot(log2(clr_summ_ds1$Mean), log2(clr_summ_ds2$Mean[m]), pch=20, xlim=c(-15,5), ylim=c(-15,5), xlab="log2(Mean CLR) DS1", ylab="log2(Mean CLR) DS2", col="slategray3", main=" ")
write.csv(cbind(clr_summ_ds1$Coordinate[use], clr_summ_ds1$Mean[use], clr_summ_ds2$Mean[m[use]]) , "SourceData/Fig5D.csv")

dev.off()
########
clr_ds1=data_DS1$circData_DS1$circCpm_filter[, -infocolsC]
Front=clr_ds1[, grep("ba9", samples1$RegionID)]
Temp=clr_ds1[, grep("ba41", samples1$RegionID)]
CB=clr_ds1[, grep("vermis", samples1$RegionID)]

Fid=samples1$BrainID[match(colnames(Front), rownames(samples1))]
Tid=samples1$BrainID[match(colnames(Temp), rownames(samples1))]
CBid=samples1$BrainID[match(colnames(CB), rownames(samples1))]


Front_matched=Front[, which(Fid%in%CBid)]
Temp_matched=Temp[, which(Tid%in%CBid)]
Fid=samples1$BrainID[match(colnames(Front_matched), rownames(samples1))]
Tid=samples1$BrainID[match(colnames(Temp_matched), rownames(samples1))]
F_index=match(Fid, CBid); T_index=match(Tid, CBid)

Cor_FCB=rep(NA, ncol(Front_matched))
for (j in c(1: length(Cor_FCB))) Cor_FCB[j]= cor(Front_matched[, j], CB[, F_index[j]], method="s", use="pair")

Cor_TCB=rep(NA, ncol(Temp_matched))
for (j in c(1: length(Cor_TCB))) Cor_TCB[j]= cor(Temp_matched[, j], CB[, T_index[j]], method="s", use="pair")

corFront=cor(Front  , method="s", use="pair")
corTemp=cor(Temp  , method="s", use="pair")
corCB=cor(CB, method="s", use="pair")
pdf("Analysis/SuppFigs/SuppFig.InterIndividual_vs_Region.pdf", height=5, width=5)
plotdata=list(F=unique(corFront),
              T=unique(corTemp),
              CB=unique(corCB),
              F_CB=Cor_FCB,
              T_CB=Cor_TCB)
boxplot(plotdata ,ylim=c(0,0.7),family="serif",
        ylab="Spearman Rho",
        col="indianred")
plotdata=list(CTX=c(unique(corFront),unique(corTemp)),
              CB=unique(corCB),
              CTX_CB=c(Cor_FCB,Cor_TCB))

boxplot(plotdata ,ylim=c(0,0.7),family="serif",
        ylab="Spearman Rho",
        col="indianred")
wilcox.test(plotdata$CTX, plotdata$CTX_CB)
wilcox.test(plotdata$CB, plotdata$CTX_CB)
wilcox.test(plotdata$CB, plotdata$CTX)
dev.off()
  
  
  
  











