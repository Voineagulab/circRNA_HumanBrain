setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
path="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/"
rm(list=ls())
library(data.table)
library(ggplot2);library(ggpubr)
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
par(pch=20, col="grey", family="serif")
############################################ LOAD circRNA DATA ############################################
exons=read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GTF/GTF_to_BED/exonsUnique.bed", sep="\t", header=FALSE)
introns=read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GTF/GTF_to_BED/introns.bed", sep="\t", header=FALSE)
#load("Analysis/DE/DE_Analysis.Rdata")
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
#Extract filtered data for genes, circRNAs and splice junctions,  DS1
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
# circCts_ds1=data_DS1$circData_DS1$circCounts_filter; 
circCpm_ds1=data_DS1$circData_DS1$circCpm_filter
# geneCts_ds1=data_DS1$geneData_DS1$counts_filter; 
geneRPKM_ds1=data_DS1$geneData_DS1$rpkm_filter
# sjMCpm_ds1=data_DS1$sjData_DS1$sjmaxCpm[which(rownames(data_DS1$sjData_DS1$sjmaxCpm)%in%rownames(circCpm_ds1)),]
# sjMCts_ds1=data_DS1$sjData_DS1$sjmaxCounts[which(rownames(data_DS1$sjData_DS1$sjmaxCounts)%in%rownames(circCpm_ds1)),]
# sjSCpm_ds1=data_DS1$sjData_DS1$sjsumCpm[which(rownames(data_DS1$sjData_DS1$sjsumCpm)%in%rownames(circCpm_ds1)),]
# sjSCts_ds1=data_DS1$sjData_DS1$sjsumCounts[which(rownames(data_DS1$sjData_DS1$sjsumCounts)%in%rownames(circCpm_ds1)),]
#Extract filtered data for genes, circRNAs and splice junctions,  DS2
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
# circCts_ds2=data_DS2$circData_DS2$circCounts_filter; 
circCpm_ds2=data_DS2$circData_DS2$circCpm_filter
# geneCts_ds2=data_DS2$geneData_DS2$counts_filter; 
geneRPKM_ds2=data_DS2$geneData_DS2$rpkm_filter
# sjMCpm_ds2=data_DS2$sjData_DS2$sjmaxCpm[which(rownames(data_DS2$sjData_DS2$sjmaxCpm)%in%rownames(circCpm_ds2)),]
# sjMCts_ds2=data_DS2$sjData_DS2$sjmaxCounts[which(rownames(data_DS2$sjData_DS2$sjmaxCounts)%in%rownames(circCpm_ds2)),]
# sjSCpm_ds2=data_DS2$sjData_DS2$sjsumCpm[which(rownames(data_DS2$sjData_DS2$sjsumCpm)%in%rownames(circCpm_ds2)),]
# sjSCts_ds2=data_DS2$sjData_DS2$sjsumCounts[which(rownames(data_DS2$sjData_DS2$sjsumCounts)%in%rownames(circCpm_ds2)),]
# load circ index data
load("DATA_TABLES/ratios_ds1.rda"); load("DATA_TABLES/ratios_ds2.rda")
ci_ds1=ratios_ds1$ci; ci_ds2=ratios_ds2$ci
# Define infocols
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)
# Select control samples 
con=which(colnames(ci_ds1)%in%samples1$Sample[which(samples1$ASD.CTL=="CTL")])
# Divide the DS1 data in CTX and CB datasets for comparison with teh AS data
circ.CBL=ci_ds1[, c(infocolsC, intersect(con, grep("vermis", colnames(ci_ds1))))]
circ.CTX=ci_ds1[, c(infocolsC, intersect(con, grep("ba", colnames(ci_ds1))))]
############################################ LOAD AND FORMAT ALT. SPLICING DATA ############################################
formatMATS=function(colData, infoData, repNames)
{
  data=data.frame(transpose(strsplit(as.character(colData), split=",")))
  colnames(data)=repNames
  data=as.data.frame(cbind(infoData, data))
  for (j in c((ncol(infoData)+1):ncol(data))) data[,j]=as.numeric(as.character(data[,j]))
  return(data)
}
names=read.csv("../DATA_TABLES/rMATS_OUTPUT/ folder names and associated input files.csv")
s1=as.character(names$SampleName[which(names$rMATS_Sample%in%"SAMPLE_1")])
s2=as.character(names$SampleName[which(names$rMATS_Sample%in%"SAMPLE_2")])
#All alternative splicing data
se=read.table("../DATA_TABLES/rMATS_OUTPUT/MATS_output/SE.MATS.JunctionCountOnly.txt", sep="\t", header=TRUE)
# a3=read.table("/Volumes/Data1/PROJECTS/circRNAs/DHG_DATA/rMATS_OUTPUT/MATS_output/A3SS.MATS.JunctionCountOnly.txt", sep="\t", header=TRUE)
# a5=read.table("/Volumes/Data1/PROJECTS/circRNAs/DHG_DATA/rMATS_OUTPUT/MATS_output/A5SS.MATS.JunctionCountOnly.txt", sep="\t", header=TRUE)
# mx=read.table("/Volumes/Data1/PROJECTS/circRNAs/DHG_DATA/rMATS_OUTPUT/MATS_output/MXE.MATS.JunctionCountOnly.txt", sep="\t", header=TRUE)
# 
# nrow(se)/(nrow(se)+nrow(a3)+nrow(a5)+nrow(mx))
# #[1] 0.7171805
# nrow(a3)/(nrow(se)+nrow(a3)+nrow(a5)+nrow(mx))
# #[1] 0.03788415
# nrow(a5)/(nrow(se)+nrow(a3)+nrow(a5)+nrow(mx))
# #[1] 0.02599131
# nrow(mx)/(nrow(se)+nrow(a3)+nrow(a5)+nrow(mx))
# #[1] 0.218944
# Based on the above, focus on Skipped Exons, i.e. cassette events
infocolsMATS=c(1:11,17:20,23);colnames(se)[infocolsMATS]
seData=list( rMATS=se,
             IJC_SAMPLE_1=formatMATS(se$IJC_SAMPLE_1, se[, infocolsMATS], s1), 
             IJC_SAMPLE_2=formatMATS(se$IJC_SAMPLE_2, se[, infocolsMATS], s2),
             
             SJC_SAMPLE_1=formatMATS(se$SJC_SAMPLE_1, se[, infocolsMATS], s1),
             SJC_SAMPLE_2=formatMATS(se$SJC_SAMPLE_2, se[, infocolsMATS], s2),
             
             IncLevel_SAMPLE_1=formatMATS(se$IncLevel1, se[, infocolsMATS], s1),
             IncLevel_SAMPLE_2=formatMATS(se$IncLevel2, se[, infocolsMATS], s2)
            )
############################################ AS data filtering ############################################
####Filter AS data for a minimum of 2 rads per inclusion junction and skipping juction respectively, in a minimum of 2 samples.
infocolsAS=c(1:16)
#CTX
keepI=which(rowSums(seData$IJC_SAMPLE_1[, -infocolsAS] >=2 ) >=2 )
keepS=which(rowSums(seData$SJC_SAMPLE_1[, -infocolsAS] >=2 ) >=2 )
keep=intersect(seData$IJC_SAMPLE_1$ID[keepI] , seData$SJC_SAMPLE_1$ID[keepS])
length(keepI);length(keepS);length(keep)

seData$IJC_SAMPLE_1=seData$IJC_SAMPLE_1[which(seData$IJC_SAMPLE_1$ID%in%keep) , ]
seData$SJC_SAMPLE_1=seData$SJC_SAMPLE_1[which(seData$SJC_SAMPLE_1$ID%in%keep) , ]
seData$IncLevel_SAMPLE_1=seData$IncLevel_SAMPLE_1[which(seData$IncLevel_SAMPLE_1$ID%in%keep) , ]
#CB
keepI=which(rowSums(seData$IJC_SAMPLE_2[, -infocolsAS] >=2 ) >=2 )
keepS=which(rowSums(seData$SJC_SAMPLE_2[, -infocolsAS] >=2 ) >=2 )
keep=intersect(seData$IJC_SAMPLE_2$ID[keepI] , seData$SJC_SAMPLE_2$ID[keepS])
seData$IJC_SAMPLE_2=seData$IJC_SAMPLE_2[which(seData$IJC_SAMPLE_2$ID%in%keep) , ]
seData$SJC_SAMPLE_2=seData$SJC_SAMPLE_2[which(seData$SJC_SAMPLE_2$ID%in%keep) , ]
seData$IncLevel_SAMPLE_2=seData$IncLevel_SAMPLE_2[which(seData$IncLevel_SAMPLE_2$ID%in%keep) , ]

pdf("Analysis/Fig4/Fig4.SplicingPlots.pdf", height=11, width=8)
par(mfrow=c(3,2))
#################################################################CTX ############################################
##format AS dataset for comparison with Circ data
se.CTX=seData$IncLevel_SAMPLE_1; colnames(se.CTX)[-infocolsAS]=formatLabels(colnames(se.CTX)[-infocolsAS])
commonSamples=intersect(colnames(se.CTX)[-infocolsAS], colnames(circ.CTX)[-infocolsC])
se.CTX$chr=paste0("chr", se.CTX$chr)
se.CTX=cbind(se.CTX[,infocolsAS], se.CTX[,commonSamples]); 
mean.f=function(x) mean(x, na.rm=TRUE)
se.CTX$mean=apply(se.CTX[,-infocolsAS], 1,mean.f); 
se.CTX$startCoord=paste(se.CTX$chr, se.CTX$exonStart_0base +1, se.CTX$strand, sep="_")
se.CTX$endCoord=paste(se.CTX$chr, se.CTX$exonEnd, se.CTX$strand, sep="_")
##format Circ dataset 
circ.CTX=cbind(circ.CTX[,infocolsC], circ.CTX[,commonSamples])
circ.CTX$mean=apply(circ.CTX[, -infocolsC], 1,mean.f)
circ.CTX$startCoord=paste(circ.CTX$Chr, circ.CTX$Start, circ.CTX$Strand, sep="_")
circ.CTX$endCoord=paste(circ.CTX$Chr, circ.CTX$End, circ.CTX$Strand, sep="_")

circ_se=which((circ.CTX$startCoord%in%se.CTX$startCoord)|(circ.CTX$endCoord%in%se.CTX$endCoord))
as_circ=which((se.CTX$startCoord%in%circ.CTX$startCoord)|(se.CTX$endCoord%in%circ.CTX$endCoord))
##########What is the proportion of skipped exons among all exons vs. circ exons
circGenes.ctx=unique(circ.CTX$EnsID)
geneRPKM_ds1$meanCTX=apply(geneRPKM_ds1[, grep("ba", colnames(geneRPKM_ds1))], 1, mean)
geneRPKM_ds1$circGenes.ctx=rownames(geneRPKM_ds1)%in%circGenes.ctx

##
exons$startCoord=paste(paste0("chr", exons[,1]), exons[,2], exons[,6], sep="_")
exons$endCoord=paste(paste0("chr", exons[,1]), exons[,3], exons[,6], sep="_")
##
exons_use=exons[which(exons[,5]%in%geneRPKM_ds1$EnsID),]
circ.CTX_use=circ.CTX[which(circ.CTX$EnsID%in%geneRPKM_ds1$EnsID),]
se.CTX_use=se.CTX[which(se.CTX$GeneID%in%geneRPKM_ds1$EnsID),]

exons_use$Circ=ifelse((exons_use$startCoord%in%circ.CTX_use$startCoord | exons_use$endCoord%in%circ.CTX_use$endCoord ), TRUE, FALSE)
exons_use$AS=ifelse((exons_use$startCoord%in%se.CTX_use$startCoord| exons_use$endCoord%in%se.CTX_use$endCoord ), TRUE, FALSE)

nExons.CTX=as.matrix(table(exons_use$Circ, exons_use$AS))
rownames(nExons.CTX)=c("non_circ", "circ"); colnames(nExons.CTX)=c("non_AS" , "AS")
write.csv(nExons.CTX, "Analysis/Fig4/Fig4.circ_skippedExonsCTX.csv")
# randomly select an equal number of genes as circ forming genes, with the same range of gene expression levels.
# From those genes do 1000 random samplings of introns with the same range of length as circ-forming introns.
# 
# b=seq(from=0, to=max(geneRPKM_ds1$meanCTX), by=0.01)
# geneRPKM_ds1$meanCTX.Bin=cut(geneRPKM_ds1$meanCTX, b)
# c=geneRPKM_ds1[which(rownames(geneRPKM_ds1)%in%circGenes.ctx), ]
# nc=geneRPKM_ds1[-which(rownames(geneRPKM_ds1)%in%circGenes.ctx), ]
# 
# nc_use=nc[sample(1:nrow(nc)) , ]
# m=match(c$meanCTX.Bin, nc_use$meanCTX.Bin)
# m=unique(m[-which(is.na(m)==TRUE)])
# result=nc_use$meanCTX[m]
# nc_use=nc_use[-m,]
# repeat {
#   m=match(c$meanCTX.Bin, nc_use$meanCTX.Bin)
#   m=unique(m[-which(is.na(m)==TRUE)])
#   result=c(result,nc_use$meanCTX[m])
#   nc_use=nc_use[-m,]  
#   print(length(m))
#   if (length(result) > nrow(c) ) {break}
# }
# result=result[c(1:nrow(c))]
# # test that you obtained a sample with similar expression distribution as the circ genes
# #par(mfrow=c(1,2))
# #boxplot(list(circ=log10(c$meanCTX) , sample=log10(result)), ylim=c(-2,4))
# #plot(density(log10(result))); lines(density(log10(c$meanCTX)))
# use_nc_genes=nc_use$EnsID[which(nc$meanCTX%in%result)]
# exons_use=exons_use[which(exons_use[,5]%in%c(use_nc_genes,circGenes.ctx)), ]
# 
# # format intron file
# introns_use=introns[which(introns$V5%in%exons_use$V5) , ] 
# colnames(introns_use)=c("Chr", "Start", "End", "Id", "Score", "Strand")
# introns_use$Chr=paste0("chr", introns_use$Chr)
# introns_use$Length=introns_use$End-introns_use$Start
# introns_use$Start=introns_use$Start-1
# introns_use$End=introns_use$End+1
# introns_use$startCoord=paste(introns_use$Chr, introns_use$Start, introns_use$Strand, sep="_")
# introns_use$endCoord=paste(introns_use$Chr, introns_use$End, introns_use$Strand, sep="_")
# introns_use_sorted=introns_use[order(introns_use$Length, decreasing=TRUE) , ]
# # compare flanking intron length for cicr_exons, as_exons and intersect.
# exons_use$leftIntronLength=introns_use_sorted$Length[match(exons_use$startCoord, introns_use_sorted$endCoord)]
# exons_use$rightIntronLength=introns_use_sorted$Length[match(exons_use$endCoord, introns_use_sorted$startCoord)]
# exons_use$leftIntronLength[which(is.na(exons_use$leftIntronLength)==TRUE)]=0
# exons_use$rightIntronLength[which(is.na(exons_use$rightIntronLength)==TRUE)]=0
# exons_use$flankingIntronLength=exons_use$leftIntronLength+exons_use$rightIntronLength
# 
# #par(mfrow=c(1,4))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$Circ, names=c("NonCirc", "Circ"), 
#         ylab="Flanking Intron Length",
#         col=c("grey", "skyblue"))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$AS, names=c("NonAS", "AS"), 
#         ylab="Flanking Intron Length",
#         col=c("grey", "indianred2"))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$AS + exons_use$Circ, 
#         names=c("Background", "AS", "Circ", "AS&Circ"), 
#         ylab="Flanking Intron Length",
#         col=c("grey", "indianred2", "skyblue", "orchid"))
# 
# #random sampling of non-circ exons with similar flanx=king intron length as circ-exons
# b=seq(from=min(exons_use$flankingIntronLength), to=max(exons_use$flankingIntronLength), by=100)
# exons_use$flankingIntronLength.Bin=cut(exons_use$flankingIntronLength, b)
# c=exons_use[which(exons_use$Circ==TRUE), ]
# nc=exons_use[which(exons_use$Circ==FALSE), ]
# ratio=array(NA, 10000)
# for(i in c(1:10000))
# {
#   s=match(c$flankingIntronLength.Bin, nc$flankingIntronLength.Bin)
#   print(i)
#   print(length(which(is.na(s)==TRUE)))
#   s=s[which(is.na(s)==FALSE)]
#   t=table(nc$AS[s])
#   print("----")
#   ratio[i]=t[2]/(t[2]+t[1])
#   nc=nc[sample(c(1:nrow(nc)), nrow(nc)),]
# }
# t=table(c$AS)
# score=t[2]/(t[2]+t[1])
# pval=length(which(ratio>score))/length(ratio)
# hist(ratio, col="lightgrey", xlim=c(0,0.3), xlab="Percent AS, CTX", sub=paste("pval=", pval))
# abline(v=score, col="red")
########## Plot distribution of Inclusion Levels for circ Exons and non-circ Exons that undergo AS
plot(lwd=3, main="AS Inc Level CTX", ylim=c(1,17),density(se.CTX$mean[-as_circ], na.rm=TRUE), col="steelblue")
lines(col="indianred",lwd=2, ylim=c(0,30), density(se.CTX$mean[as_circ], na.rm=TRUE))
wilcox.test(se.CTX$mean[-as_circ],se.CTX$mean[as_circ] )
t.test(se.CTX$mean[-as_circ],se.CTX$mean[as_circ] )
plot(0,1)
################################################################# CBL #################################################################
##format AS dataset for comparison with Circ data
se.CBL=seData$IncLevel_SAMPLE_2; colnames(se.CBL)[-infocolsAS]=formatLabels(colnames(se.CBL)[-infocolsAS])
commonSamples=intersect(colnames(se.CBL)[-infocolsAS], colnames(circ.CBL)[-infocolsC])
se.CBL$chr=paste0("chr", se.CBL$chr)
se.CBL=cbind(se.CBL[,infocolsAS], se.CBL[,commonSamples]); 

mean.f=function(x) mean(x, na.rm=TRUE)
se.CBL$mean=apply(se.CBL[,-infocolsAS], 1,mean.f); 
se.CBL$startCoord=paste(se.CBL$chr, se.CBL$exonStart_0base +1, se.CBL$strand, sep="_")
se.CBL$endCoord=paste(se.CBL$chr, se.CBL$exonEnd, se.CBL$strand, sep="_")
##format Circ dataset 
circ.CBL=cbind(circ.CBL[,infocolsC], circ.CBL[,commonSamples])
circ.CBL$mean=apply(circ.CBL[, -infocolsC], 1,mean.f)
circ.CBL$startCoord=paste(circ.CBL$Chr, circ.CBL$Start, circ.CBL$Strand, sep="_")
circ.CBL$endCoord=paste(circ.CBL$Chr, circ.CBL$End, circ.CBL$Strand, sep="_")

circ_se=which((circ.CBL$startCoord%in%se.CBL$startCoord)|(circ.CBL$endCoord%in%se.CBL$endCoord))
as_circ=which((se.CBL$startCoord%in%circ.CBL$startCoord)|(se.CBL$endCoord%in%circ.CBL$endCoord))
##########What is the proportion of skipped exons among all exons vs. circ exons
circGenes.CBL=unique(circ.CBL$EnsID)
geneRPKM_ds1$meanCBL=apply(geneRPKM_ds1[, grep("verm", colnames(geneRPKM_ds1))], 1, mean)
geneRPKM_ds1$circGenes.CBL=rownames(geneRPKM_ds1)%in%circGenes.CBL

##
exons_use=exons[which(exons[,5]%in%geneRPKM_ds1$EnsID),]
circ.CBL_use=circ.CBL[which(circ.CBL$EnsID%in%geneRPKM_ds1$EnsID),]
se.CBL_use=se.CBL[which(se.CBL$GeneID%in%geneRPKM_ds1$EnsID),]

exons_use$Circ=ifelse((exons_use$startCoord%in%circ.CBL_use$startCoord | exons_use$endCoord%in%circ.CBL_use$endCoord ), TRUE, FALSE)
exons_use$AS=ifelse((exons_use$startCoord%in%se.CBL_use$startCoord| exons_use$endCoord%in%se.CBL_use$endCoord ), TRUE, FALSE)

nExons.CBL=as.matrix(table(exons_use$Circ, exons_use$AS))
rownames(nExons.CBL)=c("non_circ", "circ"); colnames(nExons.CBL)=c("non_AS" , "AS")
write.csv(nExons.CBL, "Analysis/Fig4/Fig4.circ_skippedExonsCBL.csv")
# # randomly select an equal number of genes as circ forming genes, with the same range of gene expression levels.
# # From those genes do 1000 random samplings of introns with the same range of length as circ-forming introns.
# # 
# b=seq(from=0, to=max(geneRPKM_ds1$meanCBL), by=0.01)
# geneRPKM_ds1$meanCBL.Bin=cut(geneRPKM_ds1$meanCBL, b)
# c=geneRPKM_ds1[which(rownames(geneRPKM_ds1)%in%circGenes.CBL), ]
# nc=geneRPKM_ds1[-which(rownames(geneRPKM_ds1)%in%circGenes.CBL), ]
# 
# nc_use=nc[sample(1:nrow(nc)) , ]
# m=match(c$meanCBL.Bin, nc_use$meanCBL.Bin)
# m=unique(m[-which(is.na(m)==TRUE)])
# result=nc_use$meanCBL[m]
# nc_use=nc_use[-m,]
# repeat {
#   m=match(c$meanCBL.Bin, nc_use$meanCBL.Bin)
#   m=unique(m[-which(is.na(m)==TRUE)])
#   result=c(result,nc_use$meanCBL[m])
#   nc_use=nc_use[-m,]  
#   print(length(m))
#   if (length(result) > nrow(c) ) {break}
# }
# result=result[c(1:nrow(c))]
# # test that you obtained a sample with similar expression distribution as the circ genes
# # par(mfrow=c(1,2))
# # boxplot(list(circ=log10(c$meanCBL) , sample=log10(result)), ylim=c(-2,4))
# # plot(density(log10(result))); lines(density(log10(c$meanCBL)))
# use_nc_genes=nc_use$EnsID[which(nc$meanCBL%in%result)]
# exons_use=exons_use[which(exons_use[,5]%in%c(use_nc_genes,circGenes.CBL)), ]
# # format intron file
# introns_use=introns[which(introns$V5%in%exons_use$V5) , ] 
# colnames(introns_use)=c("Chr", "Start", "End", "Id", "Score", "Strand")
# introns_use$Chr=paste0("chr", introns_use$Chr)
# introns_use$Length=introns_use$End-introns_use$Start
# introns_use$Start=introns_use$Start-1
# introns_use$End=introns_use$End+1
# introns_use$startCoord=paste(introns_use$Chr, introns_use$Start, introns_use$Strand, sep="_")
# introns_use$endCoord=paste(introns_use$Chr, introns_use$End, introns_use$Strand, sep="_")
# introns_use_sorted=introns_use[order(introns_use$Length, decreasing=TRUE) , ]
# # compare flanking intron length for cicr_exons, as_exons and intersect.
# exons_use$leftIntronLength=introns_use_sorted$Length[match(exons_use$startCoord, introns_use_sorted$endCoord)]
# exons_use$rightIntronLength=introns_use_sorted$Length[match(exons_use$endCoord, introns_use_sorted$startCoord)]
# exons_use$leftIntronLength[which(is.na(exons_use$leftIntronLength)==TRUE)]=0
# exons_use$rightIntronLength[which(is.na(exons_use$rightIntronLength)==TRUE)]=0
# exons_use$flankingIntronLength=exons_use$leftIntronLength+exons_use$rightIntronLength
# 
# #par(mfrow=c(1,4))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$Circ, names=c("NonCirc", "Circ"), main="CB",
#         ylab="Flanking Intron Length",
#         col=c("grey", "skyblue"))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$AS, names=c("NonAS", "AS"), main="CB",
#         ylab="Flanking Intron Length",
#         col=c("grey", "indianred2"))
# boxplot(log10(exons_use$flankingIntronLength) ~ exons_use$AS + exons_use$Circ, main="CB",
#         names=c("Background", "AS", "Circ", "AS&Circ"), 
#         ylab="Flanking Intron Length",
#         col=c("grey", "indianred2", "skyblue", "orchid"))
# 
# #random sampling of non-circ exons with similar flanx=king intron length as circ-exons
# b=seq(from=min(exons_use$flankingIntronLength), to=max(exons_use$flankingIntronLength), by=100)
# exons_use$flankingIntronLength.Bin=cut(exons_use$flankingIntronLength, b)
# c=exons_use[which(exons_use$Circ==TRUE), ]
# nc=exons_use[which(exons_use$Circ==FALSE), ]
# ratio=array(NA, 10000)
# for(i in c(1:10000))
# {
#   s=match(c$flankingIntronLength.Bin, nc$flankingIntronLength.Bin)
#   print(i)
#   print(length(which(is.na(s)==TRUE)))
#   s=s[which(is.na(s)==FALSE)]
#   t=table(nc$AS[s])
#   print("----")
#   ratio[i]=t[2]/(t[2]+t[1])
#   nc=nc[sample(c(1:nrow(nc)), nrow(nc)),]
# }
# 
# t=table(c$AS)
# score=t[2]/(t[2]+t[1])
# pval=length(which(ratio>score))/length(ratio)
# hist(ratio, col="lightgrey", xlim=c(0,0.3), xlab="Percent AS, CB", sub=paste("pval=", pval))
# abline(v=score, col="red")
########## Plot distribution of Inclusion Levels for circ Exons and non-circ Exons that undergo AS
plot(lwd=3, main="AS Inc Level CBL", ylim=c(1,17),density(se.CBL$mean[-as_circ], na.rm=TRUE), col="steelblue")
lines(col="indianred",lwd=2, ylim=c(0,30), density(se.CBL$mean[as_circ], na.rm=TRUE))
wilcox.test(se.CBL$mean[-as_circ],se.CBL$mean[as_circ] )
t.test(se.CBL$mean[-as_circ],se.CBL$mean[as_circ] )
dev.off()
######################################################### N circ isoforms vs N AS exons #################################################################
total_nExons = as.data.frame(table(exons[, 5]))
#N circ RNA isoforms
x=circCpm_ds1$Symbol[which(circCpm_ds1$Id %in% circCpm_ds2$Id)]
hs1 = as.data.frame(table(x))
colnames(hs1) = c("Symbol", "nCirc")
#hs1 = hs1[-which(hs1$Symbol %in% "NA"), ]
hs1 = hs1[-which(hs1$nCirc ==0 ), ]
table(hs1$nCirc >= 5)
#FALSE  TRUE
#3025   446

#N AS exons
as=as.data.frame(table(seData$IncLevel_SAMPLE_1$geneSymbol))

#Combine in the same dataframe
hs1$nASExons=as$Freq[match(hs1$Symbol, as$Var1)]

#Replace NA entries for AS exons with 0
hs1$nASExons[which(is.na(hs1$nASExons)==TRUE)]=0

#Plot
pdf("Analysis/Fig4/Fig4.Hotspotgenes_vs_AS.pdf", height=6, width=8)
plot(hs1$nCirc[which(hs1$nCirc>=10)], hs1$nASExons[which(hs1$nCirc>=10)], pch=20, xlab="Number of circRNAs per gene", ylab="Number of AS events per gene",
     sub=paste("rho=", round(cor(hs1$nCirc ,hs1$nASExons , method="s"),2)),
     xlim=c(10,50))
text(pos=3,cex=0.7,hs1$nCirc[which(hs1$nCirc>=15)], hs1$nASExons[which(hs1$nCirc>=15)], hs1$Symbol[which(hs1$nCirc>=15)])
dev.off()
# Plot in large size for gene labels
pdf("Analysis/Fig4/Fig4.Hotspotgenes_vs_AS_large.pdf", height=15, width=20)
plot(hs1$nCirc[which(hs1$nCirc>=10)], hs1$nASExons[which(hs1$nCirc>=10)], pch=20, xlab="Number of circRNAs per gene", ylab="Number of AS events per gene",
     sub=paste("rho=", round(cor(hs1$nCirc ,hs1$nASExons , method="s"),2)),
     xlim=c(10,50))
text(pos=3,cex=0.7,hs1$nCirc[which(hs1$nCirc>=15)], hs1$nASExons[which(hs1$nCirc>=15)], hs1$Symbol[which(hs1$nCirc>=15)])
dev.off()

# Add and correct for total N exons
nExons_Gene=table(exons$V5)
hs1$EnsID=circCpm_ds1$EnsID[match(hs1$Symbol, circCpm_ds1$Symbol)]
hs1$Total_nExons=nExons_Gene[match(hs1$EnsID, names(nExons_Gene))]
hs1$nCirc_prop=as.numeric(as.character(hs1$nCirc/hs1$Total_nExon))
hs1$nAS_prop=as.numeric(as.character(hs1$nASExons/hs1$Total_nExon))
#hs1=hs1[-which(is.na(hs1$nAS_prop)==TRUE) ,]
# Cor test
cor.test(hs1$nCirc , hs1$nASExons , method="s")

elim=which(hs1$Symbol%in%c("PTK2", "ATRNL1"))
cor.test(hs1$nCirc[-elim] , hs1$nASExons[-elim] , method="s")

# Spearman's rank correlation rho
# 
# data:  hs1$nCirc and hs1$nASExons
# S = 4.493e+09, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.3225178 

cor.test(hs1$nCirc_prop , hs1$nAS_prop , method="s")
# Spearman's rank correlation rho
# 
# data:  hs1$nCirc_prop and hs1$nAS_prop
# S = 6019300000, p-value = 7.107e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.1363629 