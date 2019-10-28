rm(list=ls())
library(WGCNA)
library(ggplot2);
library(ggpubr); 
library(dplyr);
library(data.table)
library(goseq)
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

############################################ CALCULATE SUMMARY INFO ACROSS DATASETS #################
#genes rpkm
geneRPKM_summ_ds1=calcSummaryData(geneRPKM_ds1, infocolsG, samples1,th=1)
geneRPKM_summCTL_ds1=calcSummaryData(geneRPKM_ds1, infocolsG, samples1[which(samples1$ASD.CTL=="CTL"),],th=1)
geneRPKM_summ_ds2=calcSummaryData(geneRPKM_ds2,infocolsG, samples2,th=1)
geneRPKM_summCTL_ds2=calcSummaryData(geneRPKM_ds2,infocolsG, samples2[which(samples2$ASD.CTL=="CTL"),],th=1)

#circ
circCts_summ_ds1=calcSummaryData(circCts_ds1, infocolsC, samples1, th=2)
circCts_summCTL_ds1=calcSummaryData(circCts_ds1, infocolsC, samples1[which(samples1$ASD.CTL=="CTL"),], th=2)
circCpm_summ_ds1=calcSummaryData(circCpm_ds1, infocolsC, samples1, th=0.1)
circCpm_summCTL_ds1=calcSummaryData(circCpm_ds1, infocolsC, samples1[which(samples1$ASD.CTL=="CTL"),], th=0.1)

circCts_summ_ds2=calcSummaryData(circCts_ds2, infocolsC, samples2, th=2)
circCts_summCTL_ds2=calcSummaryData(circCts_ds2, infocolsC, samples2[which(samples2$ASD.CTL=="CTL"),], th=2)
circCpm_summ_ds2=calcSummaryData(circCpm_ds2, infocolsC, samples2, th=0.1)
circCpm_summCTL_ds2=calcSummaryData(circCpm_ds2, infocolsC, samples2[which(samples2$ASD.CTL=="CTL"),], th=0.1)
#sj max
sjMCpm_summ_ds1=calcSummaryData(sjMCpm_ds1, infocolsS, samples1, 0.1)
sjMCpm_summCTL_ds1=calcSummaryData(sjMCpm_ds1, infocolsS, samples1[which(samples1$ASD.CTL=="CTL"),], 0.1)
sjMCpm_summ_ds2=calcSummaryData(sjMCpm_ds2, infocolsS, samples2, 0.1)
sjMCpm_summCTL_ds2=calcSummaryData(sjMCpm_ds2, infocolsS, samples2[which(samples2$ASD.CTL=="CTL"),], 0.1)
# CI
ci_summ_ds1=calcSummaryData(ci_ds1, infocolsC, samples1, 0.1)
ci_summCTL_ds1=calcSummaryData(ci_ds1, infocolsC, samples1[which(samples1$ASD.CTL=="CTL"),], 0.01)
ci_summ_ds2=calcSummaryData(ci_ds2, infocolsC, samples2, 1)
ci_summCTL_ds2=calcSummaryData(ci_ds2, infocolsC, samples2[which(samples2$ASD.CTL=="CTL"),], 0.01)
# CLR
clr_summ_ds1=calcSummaryData(clr_ds1, infocolsC, samples1, 1)
clr_summCTL_ds1=calcSummaryData(clr_ds1, infocolsC, samples1[which(samples1$ASD.CTL=="CTL"),], 0.01)
clr_summ_ds2=calcSummaryData(clr_ds2, infocolsC, samples2, 1)
clr_summCTL_ds2=calcSummaryData(clr_ds2, infocolsC, samples2[which(samples2$ASD.CTL=="CTL"),], 0.01)
#
nGenePerSample1=colSums(geneRPKM_ds1[, -infocolsG] >= 1)
nGenePerSample2=colSums(geneRPKM_ds2[, -infocolsG] >= 1)
nCircPerSample1=colSums(circCpm_ds1[, -infocolsC] >= 0.1)
nCircPerSample2=colSums(circCpm_ds2[, -infocolsC] >= 0.1)
#nCircPerSample1=colSums(circCts_ds1[, -infocolsC] >= 2)
#nCircPerSample2=colSums(circCts_ds2[, -infocolsC] >= 2)
############################################ COMPARE WITH CIRC BASE ############################################
# Load and format human circBase
circBase=read.delim("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/PublishedDatasets/hsa_hg19_circRNA_download01.2018.txt", header = TRUE, sep="\t")
# Note: the "All huamn circRNAs" file from circBase does NOT include the data from Rybak-Wolf et al., so these had to be added:
rw=read.delim("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/PublishedDatasets/hsa_hg19_Rybak2015_download01.2018.txt", header = TRUE, sep="\t")
circBase=rbind(rw,circBase)
colnames(circBase)[1] <- c("Chr")
Coordinate=paste(circBase$Chr,circBase$start +1 , circBase$end,  circBase$strand, sep="_")
Id=paste(circBase$Chr,circBase$start +1 , circBase$end,  sep="_")
circBase <- cbind(Coordinate, Id, circBase)
circBase=circBase[match(unique(circBase$Coordinate), circBase$Coordinate),]
circBase_Brain=circBase[grep("Rybak", circBase$circRNA.study), ]
write.csv(circBase, "SourceData/circBase.csv")
#Note: because or the order of rbind followed by match, the circRNAs common between circBase and RW will be labeled as "Rybak", which is what we want so that we can trace their exp in brain
circCts_summ_ds1$circ_type="Novel";circCts_summ_ds1$circ_type[which(circCts_summ_ds1$Id%in%circBase$Id)]='circBase';circCpm_summ_ds1$circ_type="Novel";circCpm_summ_ds1$circ_type[which(circCpm_summ_ds1$Id%in%circBase$Id)]='circBase'
circCts_summ_ds2$circ_type="Novel";circCts_summ_ds2$circ_type[which(circCts_summ_ds2$Id%in%circBase$Id)]='circBase'
circCpm_summ_ds2$circ_type="Novel";circCpm_summ_ds2$circ_type[which(circCpm_summ_ds2$Id%in%circBase$Id)]='circBase'

############################################ DATASET OVERVIEW #####################################################
#################  DS1 and DS2 Sample Info 
t1=table(samples1$ASD.CTL , samples1$Region)
write.csv(t1, "Analysis/Fig1/Fig1.sampletable_DS1.csv")
t2=table(samples2$ASD.CTL , samples2$Region)
write.csv(t2, "Analysis/Fig1/Fig1.sampletable_DS2.csv")

################# DS1 and DS2 Sample Info (SuppTables)
savecols=colnames(samples1)[c(1:3, 5:12, 55:61)]
write.csv(samples1[, savecols], "Analysis/SuppTables/SuppTable_SampleInfo_DS1.csv",row.names=FALSE)
write.csv(samples2[,savecols], "Analysis/SuppTables/SuppTable_SampleInfo_DS2.csv",row.names=FALSE)

################# DS1 and DS2 circRNA data with CircBase annotation (SuppTables)
useCircBase=circBase[, c(2,7:14)]
colnames(useCircBase)[-1]=paste("circBase", colnames(useCircBase)[-1], sep="_")
SuppTable_DS1_annot=circCpm_summ_ds1; SuppTable_DS2_annot=circCpm_summ_ds2; 
SuppTable_DS1_annot=merge(y=useCircBase, x=SuppTable_DS1_annot, all.x=TRUE, by.x="Id")
SuppTable_DS2_annot=merge(y=useCircBase, x=SuppTable_DS2_annot, all.x=TRUE, by.x="Id")
rownames(SuppTable_DS1_annot)=SuppTable_DS1_annot$Coordinate;rownames(SuppTable_DS2_annot)=SuppTable_DS2_annot$Coordinate;
# keepcols=c(2:9,23,24, 26:34)
# SuppTable_DS1_annot=SuppTable_DS1_annot[,keepcols];SuppTable_DS2_annot=SuppTable_DS2_annot[,keepcols]
SuppTable_DS1_exp=data_DS1$circData_DS1$circCpm_filter[, -infocolsC];
SuppTable_DS2_exp=data_DS2$circData_DS2$circCpm_filter[, -infocolsC]
write.csv(SuppTable_DS1_annot, "Analysis/SuppTables/SuppTable_DS1_annot.csv")
write.csv(SuppTable_DS2_annot, "Analysis/SuppTables/SuppTable_DS2_annot.csv")
write.csv(SuppTable_DS1_exp, "Analysis/SuppTables/SuppTable_DS1_exp.csv")
write.csv(SuppTable_DS2_exp, "Analysis/SuppTables/SuppTable_DS2_exp.csv")

################# RIMS2 example
r1=SuppTable_DS1_annot[grep("RIMS2", SuppTable_DS1_annot$Symbol),]
r2=SuppTable_DS2_annot[grep("RIMS2", SuppTable_DS2_annot$Symbol),]
r1=r1[intersect(rownames(r1), rownames(r2)),];r2=r2[intersect(rownames(r1), rownames(r2)),]
dim(r1);dim(r2);
#15 common RIMS circRNAsbetween DS1 and DS2
table(is.na(r1$circBase_circRNA.ID))
#7 novel circRNAs
r1=SuppTable_DS1_annot[intersect(grep("RIMS2", SuppTable_DS1_annot$Symbol), which(SuppTable_DS1_annot$nexpPerFeature>=2 )),]
r2=SuppTable_DS2_annot[intersect(grep("RIMS2", SuppTable_DS2_annot$Symbol), which(SuppTable_DS2_annot$nexpPerFeature>=2 )),]
r1=r1[intersect(rownames(r1), rownames(r2)),];r2=r2[intersect(rownames(r1), rownames(r2)),]
dim(r1);dim(r2)
table(is.na(r1$circBase_circRNA.ID))
#7 highly expressed , 3 of these are novel
write.csv(r1,"Analysis/Fig2/Fig2.RIMS2_ds1.csv")
write.csv(r2,"Analysis/Fig2/Fig2.RIMS2_ds2.csv")

################# Gene Ontology Analysis of CircGenes 
# How many genes form circ RNAs?
length(unique(circCpm_summ_ds1$EnsID))
#[1] 4555
length(unique(circCpm_summ_ds2$EnsID))
#[1] 3650
circGenes1=unique(circCpm_ds1$EnsID); length(circGenes1) ; circGenes2=unique(circCpm_ds2$EnsID); length(circGenes2); length(intersect(circGenes1 , circGenes2))
expGenes1=unique(geneRPKM_ds1$EnsID); length(expGenes1); expGenes2=unique(geneRPKM_ds2$EnsID); length(expGenes2);length(intersect(expGenes1 , expGenes2))
expGenes=intersect(expGenes1, expGenes2); circGenes=intersect(circGenes1,circGenes2)
circGenes_goData=ifelse(expGenes%in%circGenes, 1,0); 
names(circGenes_goData)= expGenes

## Run gProfiler
library(gprofiler2)
go_path=gost(circGenes, organism = "hsapiens", ordered_query = FALSE,
             multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
             measure_underrepresentation = FALSE, evcodes = FALSE,
             user_threshold = 0.05, correction_method ="bonferroni", custom_bg = expGenes,
             numeric_ns = "", sources = "KEGG")
kegg=go_path$result
write.csv(kegg, "Analysis/SuppTables/SuppTable.circGenes_kegg_gProfiler.csv")
## run GoSeq
circGenes_goResults=runGoSeq(circGenes_goData); circGenes_goResults=as.data.frame(circGenes_goResults)
circGenes_goResults=circGenes_goResults[which(circGenes_goResults$FDR < 0.05) , ]
circGenes_goResults$EnsIDs=NA
circGenes_goResults$Symbols=NA

genes2go=getgo(names(circGenes_goData),'hg19','ensGene') 
go2genes=goseq:::reversemapping(genes2go)

for (j in c(1:nrow(circGenes_goResults)))
  {
  #print(j)
  t=match(circGenes_goResults$category[j], names(go2genes))
  ids=intersect(names(circGenes_goData)[circGenes_goData==1], go2genes[[t]])
  symbols=geneRPKM_ds1$Symbol[match(ids, geneRPKM_ds1$EnsID)]
  #circGenes_goResults$EnsIDs[j]=toString(ids);circGenes_goResults$EnsIDs[j]=gsub(", ", ";", circGenes_goResults$EnsIDs[j])
  circGenes_goResults$Symbols[j]=toString(symbols);circGenes_goResults$Symbols[j]=gsub(", ", ";", circGenes_goResults$Symbols[j])
}
write.csv(circGenes_goResults, "Analysis/Fig2/Fig2_and_SuppTable.circGenes_goResults.csv", row.names=FALSE)

## run additional enrichment analyses . Note: this is done using symbols 
calculate.overlaps=function(files)
{
  #each file in the list should contain a gene category. The function calculates the significance of the overlap with circRNA expressing genes.
  #hyper-geometric test p-values are Bonferroni corrected.
  p.hyper=rep(0, length(files)); prop=rep(0, length(files)); prop.exp=rep(0, length(files)); mg=rep(" ", length(files))
  names(prop)=names(prop.exp)=gsub("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/", "", files, fixed=TRUE)
  for (j in c(1:length(files)))
  {
    g=read.table(files[j], header=TRUE)
    p.hyper[j]=phyper(length(intersect(g[,1], circGenes_Symbol)),
                      length(intersect(g[,1], geneR)),
                      length(geneR)-length(intersect(g[,1], geneR)),
                      length(circGenes_Symbol),
                      lower.tail = FALSE)
    prop[j]=length(intersect(g[,1], circGenes_Symbol))/length(g[,1])
    prop.exp[j]=length(intersect(g[,1], circGenes_Symbol))/length(intersect(g[,1], geneR))
    mg[j]=toString(intersect(g[,1], circGenes_Symbol))
  }
  ov=data.frame(p.hyper*length(files), prop,prop.exp, mg); colnames(ov)=c("p.hyper.bonf", "prop", "prop.exp", "Overlap_Genes")
  return(ov)
}

##ASD, ID, SCZ genes
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/ensGeneTxSymbol.02.2019.rda")
circGenes_Symbol=ens$Symbol[match(circGenes, ens$EnsID)]
geneR=intersect(data_DS1$geneData_DS1$rpkm_filter$Symbol, data_DS2$geneData_DS2$rpkm_filter$Symbol)
labelR=ifelse(geneR%in%circGenes_Symbol, "circGenes","background"); names(labelR)=circGenes_Symbol
ndd_files=paste0("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/", c("ID_Vissers_et_al. 2016_supp_table.txt", "SFARI-15-06-2018_Syndromic_or_Score_1_4.txt", "SZDB_minScore2.txt"))
ndd_overlaps=calculate.overlaps(ndd_files)
write.csv(ndd_overlaps, "Analysis/circGenes_NDD_overlaps.csv")
##Neuronal subtype markers
ex=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/Lake_ST3_ExcitatoryNeurons.csv")
inh=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/Lake_ST3_InhibitoryNeurons.csv")

path="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/"
ex$Cluster=as.character(ex$Cluster)
exN=unique(ex$Cluster)
for (j in c(1:length(exN)))
  {
  mk=ex$Gene[which((ex$Cluster%in%exN[j])&(ex$Average.Difference..log.fold.change. > 1))]
  if (length(mk) >0) write.table(mk, sep="\t", row.names=FALSE, col.names=exN[j], file=paste0(path, "Lake_", as.character(exN[j]), "_markers.txt"), quote=FALSE)
  }
inhN=unique(as.character(inh$Cluster));
for (j in c(1:length(inhN)))
{
  mk=inh$Gene[which((inh$Cluster%in%inhN[j])&(inh$Average.Difference..log.fold.change. > 1))]
  if (length(mk) >0) write.table(mk, sep="\t", row.names=FALSE, col.names=inhN[j], file=paste0(path, "Lake_", as.character(inhN[j]), "_markers.txt"), quote=FALSE)
}
lake_files=paste0(path, "Lake_", c(exN, inhN), "_markers.txt")
lake_files=lake_files[-3] # Ex3a did not have any markes with log.FC>1

neuronal_marker_overlaps=calculate.overlaps(lake_files)
neuronal_marker_overlaps=neuronal_marker_overlaps[which((neuronal_marker_overlaps$p.hyper.bonf <0.05)&(neuronal_marker_overlaps$prop >0)), ]
write.csv(neuronal_marker_overlaps, "Analysis/SuppTables/SuppTables.circGenes_neuronalMarker_overlaps.csv")

############################################ Novel CircRNAS  #####################################################
novel1=which(circCts_summ_ds1$circ_type%in%"Novel");novel2=which(circCts_summ_ds2$circ_type%in%"Novel")
length(novel1);length(novel2); length(intersect(circCts_summ_ds1$Coordinate[novel1],circCts_summ_ds2$Coordinate[novel2] ))/length(novel2)
# [1] 1548
# [1] 692
# [1] 0.8309249
################# CircRNA Genomic Annotation tables 
annotTableDS1=table(circCts_summ_ds1$summaryAnnot); annotTableDS1N=table(circCts_summ_ds1$summaryAnnot[novel1])
annotTableDS2=table(circCts_summ_ds2$summaryAnnot); annotTableDS2N=table(circCts_summ_ds2$summaryAnnot[novel2])
annotPlot=rbind(annotTableDS1,annotTableDS1N, annotTableDS2, annotTableDS2N)
write.csv(annotPlot, "Analysis/SuppFigs/SuppFig.annotPlot.csv")

################# Plot Venn diagrams for comparison with circBase
pdf("Analysis/Fig2/Fig2.vennDiagrams_DS1_DS2_circBase.pdf", height=3, width=3)
plotVennCirc(circCpm_ds1, circCpm_ds2, circBase)
dev.off()

################# Properties of novel circRNAs 
p1 <- ggplot(circCts_summ_ds1, aes(x=circ_type, y=nexpPerFeature, fill=circ_type)) +
  geom_violin() + geom_boxplot(width=0.1, fill="white") +
  ggtitle("DS1") + ylab("N Samples Exp > 2Counts")
p2 <- ggplot(circCts_summ_ds2, aes(x=circ_type, y=nexpPerFeature, fill=circ_type)) + 
  geom_violin() + geom_boxplot(width=0.1, fill="white") +
  ggtitle("DS2") + ylab("N Samples Exp > 2Counts")
p3 <- ggplot(circCpm_summ_ds1, aes(x=circ_type, y=log2(Mean), fill=circ_type)) +
  geom_violin() + geom_boxplot(width=0.1, fill="white") +
  ggtitle("DS1") + ylab("Mean Cpm") + scale_fill_manual(values=c("steelblue1", "yellow2"))
p4 <- ggplot(circCpm_summ_ds2, aes(x=circ_type, y=log2(Mean), fill=circ_type)) +
  geom_violin() + geom_boxplot(width=0.1, fill="white") +
  ggtitle("DS2") + ylab("Mean Cpm") + scale_fill_manual(values=c("steelblue1", "yellow2"))
fig=ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, labels=c("A", "B", "C", "D"))
pdf("Analysis/SuppFigs/SuppFig.NovelCirc_Properties.pdf", height=8, width=8)
print(fig)
dev.off()
################# Properties of novel circRNAs: N Samples Exp 
# Novel circRNAs expressed in at least 10 samples
n1=which(circCts_summ_ds1$nexpPerFeature >= 10); n2=which(circCts_summ_ds2$nexpPerFeature >= 10)
pdf("Analysis/Fig2/Fig2.NovelCirc_nSamplesExp.pdf", height=8, width=8)
par(mfrow=c(3,3))
# select novel circRNAs expressed in at least 10 samples; plot hist of nsamples in which these circRNAs are expressed
hist(circCts_summ_ds1$nexpPerFeature[intersect(n1, novel1)], col="mediumorchid", breaks=20, main="", xlab="Number of samples", ylab="N circRNAs")
hist(circCts_summ_ds2$nexpPerFeature[intersect(n2, novel2)], col="seagreen3", breaks=20,  main="", xlab="Number of samples", ylab="N circRNAs")
dev.off()
write.csv(circCts_summ_ds1, "Analysis/Fig2/Fig2.circCts_summ_ds1.csv")
write.csv(circCts_summ_ds2, "Analysis/Fig2/Fig2b.circCts_summ_ds2.csv")
length(intersect(novel1,n1)); length(intersect(novel2,n2))
#[1] 727
#[1] 287
#How many come from novel circRNA producing genes?
novel_g_d1=setdiff(circCts_summ_ds1$Symbol, circBase$gene.symbol)
#length(novel_g_d1)
#[1] 342
novel_g_d2=setdiff(circCts_summ_ds2$Symbol, circBase$gene.symbol)
#length(novel_g_d2)
#[1] 217
novel_g=intersect(novel_g_d1, novel_g_d2)
#length(novel_g)
#[1] 207
novel_g=cbind(novel_g, circCpm_ds1$EnsID[match(novel_g, circCpm_ds1$Symbol)])
novel_g=as.data.frame(novel_g)
colnames(novel_g)=c("Symbol","EnsID")
length(which(circCpm_ds1$Symbol%in%novel_g$Symbol));length(which(circCpm_ds2$Symbol%in%novel_g$Symbol))
#[1] 1080
#[1] 662

# What proportion of circRNAs are produced from genes detected as expressed?
length(intersect(unique(circCpm_ds1$EnsID),unique(geneRPKM_ds1$EnsID) ))/length(unique(circCpm_ds1$EnsID))
length(intersect(unique(circCpm_ds2$EnsID),unique(geneRPKM_ds2$EnsID) ))/length(unique(circCpm_ds2$EnsID))
#[1]  0.9413831
#[1]  0.9652055
# Looks reasonable.

############## CircRNA vs Linear Gene ExpressionScatterplots 
pdf("Analysis/SuppFigs/SuppFig.Circ_vs_geneExp_CTL.pdf", height=8, width=8)
par(mfrow=c(3,3))
use_ds1=which((sjMCpm_summCTL_ds1$Mean >0)&(circCpm_summCTL_ds1$Mean>0))
scatterplot.mean.f(circCpm_summCTL_ds1[use_ds1,], sjMCpm_summCTL_ds1[use_ds1, ], xlab="Circ Junction Mean Expression", ylab="Linear Junction Mean Expression", title="DS1", matchcol="Id")
scatterplot.mean.f(circCpm_summCTL_ds1[use_ds1,], geneRPKM_summCTL_ds1, xlab="Circ Junction Mean Expression", ylab="Gene-level Mean Expression", title="DS1", matchcol="Symbol")
scatterplot.mean.f(sjMCpm_summCTL_ds1[use_ds1, ], geneRPKM_summCTL_ds1, xlab="Linear Junction Mean Expression", ylab="Gene-level Mean Expression", title="DS1", matchcol="Symbol")

use_ds2=which((sjMCpm_summCTL_ds2$Mean >0)&(circCpm_summCTL_ds2$Mean>0))
scatterplot.mean.f(circCpm_summCTL_ds2[use_ds2,], sjMCpm_summCTL_ds2[use_ds2, ], xlab="Circ Junction Mean Expression", ylab="Linear Junction Mean Expression", title="DS2", matchcol="Id")
scatterplot.mean.f(circCpm_summCTL_ds2[use_ds2,], geneRPKM_summCTL_ds2, xlab="Circ Junction Mean Expression", ylab="Gene-level Mean Expression", title="DS2", matchcol="Symbol")
scatterplot.mean.f(sjMCpm_summCTL_ds2[use_ds2, ], geneRPKM_summCTL_ds2, xlab="Linear Junction Mean Expression", ylab="Gene-level Mean Expression", title="DS2", matchcol="Symbol")
plotdataDS1=list(circGenes=log2(geneRPKM_summCTL_ds1$Mean[which(geneRPKM_summCTL_ds1$EnsID%in%circGenes1)]),
                 non_circExpGenes=log2(geneRPKM_summCTL_ds1$Mean[-which(geneRPKM_summCTL_ds1$EnsID%in%circGenes1)]))
boxplot(plotdataDS1, col=c("lightgrey", "lemonchiffon4"),  notch=TRUE, main="DS1", ylab="Gene-level Mean Expression")          

plotdataDS2=list(circGenes=log2(geneRPKM_summCTL_ds2$Mean[which(geneRPKM_summCTL_ds2$EnsID%in%circGenes2)]),
                 non_circExpGenes=log2(geneRPKM_summCTL_ds2$Mean[-which(geneRPKM_summCTL_ds2$EnsID%in%circGenes2)]))
boxplot(plotdataDS2, main="DS2",col=c("lightgrey", "lemonchiffon4"),  notch=TRUE, ylab="Gene-level Mean Expression")  

dev.off()

############## Hotspot genes
hs1=as.data.frame(table(circCpm_ds1$EnsID))
length(which((hs1$Freq >= 5)))
#[1] 878
hs2=as.data.frame(table(circCpm_ds2$EnsID))
length(which((hs2$Freq >= 5)))
#[1] 509
novel_g$hotSpot1=hs1$Freq[match(novel_g$EnsID ,hs1$Var1)]
novel_g$hotSpot2=hs2$Freq[match(novel_g$EnsID ,hs2$Var1)]
write.csv(novel_g, "Analysis/NovelCircRNAGenes_Hotspots.csv")

pdf("Analysis/Fig3/Fig3.hotspotGenesCTL.pdf", height=8, width=8)
par(mfrow=c(3,2))
par(pch=20, col="grey", family="serif")
####DS1
ncirc_per_gene1=table(circCpm_summCTL_ds1$EnsID[-grep("NA", circCpm_summCTL_ds1$EnsID)])
write.csv(ncirc_per_gene1, "Analysis/Fig3/Fig3.ncirc_per_gene1.csv")
hist(ncirc_per_gene1, col="tomato", breaks=seq(from=0 , to=60, by=1), xlab="N circRNA isoforms", ylab="N genes", main="DS1")
hist(ncirc_per_gene1[which(ncirc_per_gene1 >= 5)], col="tomato",xlim=c(5,60),  breaks=seq(from=0 , to=60, by=1), xlab="N circRNA isoforms", ylab="N genes", main="DS1")
####DS2
ncirc_per_gene2=table(circCpm_summCTL_ds2$EnsID[-grep("NA", circCpm_summCTL_ds2$EnsID)])
hist(ncirc_per_gene2, col="tomato", breaks=seq(from=1 , to=60, by=1), xlab="N circRNA isoforms", ylab="N genes", main="DS2")
hist(ncirc_per_gene2[which(ncirc_per_gene2 >= 5)], col="tomato", xlim=c(5,60), breaks=seq(from=1 , to=60, by=1), xlab="N circRNA isoforms", ylab="N genes", main="DS2")
dev.off()

############## Major Isoform Expression
#pdf("Analysis/DatasetCharacterisationPlots/Fig3.majorIsoform_CTL.pdf", height=9, width=7.5)
#par(mfrow=c(2,1));par(pch=20, col="grey", family="serif")
####DS1
#CPM
# Calculate the Max and Sum of expression for each circRNA
na1=grep("NA", circCpm_summCTL_ds1$EnsID)
circ_ds1_agMax=aggregate(circCpm_summCTL_ds1$Mean[-na1] ~ circCpm_summCTL_ds1$EnsID[-na1], FUN = max)
circ_ds1_agSum=aggregate(circCpm_summCTL_ds1$Mean[-na1] ~ circCpm_summCTL_ds1$EnsID[-na1], FUN = sum)
# Match the order of circRNAs in the expression dataframes, and the ncirc_per_gene vector
circ_ds1_agMax=circ_ds1_agMax[match(names(ncirc_per_gene1), circ_ds1_agMax[,1]), ]
circ_ds1_agSum=circ_ds1_agSum[match(names(ncirc_per_gene1), circ_ds1_agSum[,1]), ]
# Calculate major isoform relative expression 
major_isoform_ratio=circ_ds1_agMax[,2]/circ_ds1_agSum[,2]
# Boxplot
#boxplot(major_isoform_ratio ~ ncirc_per_gene1, main="Circ Cpm DS1", ylim=c(0,1), col="lightgrey", ,xlab="N circRNA isoforms per gene", ylab="Major Isoform Relative Expression")

#x=c(1:length(table(ncirc_per_gene1))); y=1/x; lines(x,y, col="red")
#write.csv(cbind(major_isoform_ratio, ncirc_per_gene1), "SourceData/Fig3C.MajorIsoform.csv")
#CLR
# Calculate the Max and Sum of CLR for each circRNA
na1=grep("NA", clr_summCTL_ds1$EnsID)
clr_ds1_agMax=aggregate(clr_summCTL_ds1$Mean[-na1] ~ clr_summCTL_ds1$EnsID[-na1], FUN = max)
clr_ds1_agSum=aggregate(clr_summCTL_ds1$Mean[-na1] ~ clr_summCTL_ds1$EnsID[-na1], FUN = sum)
# Match the order of circRNAs in the CLR dataframes, and the nclr_per_gene vector
nclr_per_gene=table(clr_summCTL_ds1$EnsID[-na1])
clr_ds1_agMax=clr_ds1_agMax[match(names(nclr_per_gene), clr_ds1_agMax[,1]), ]
clr_ds1_agSum=clr_ds1_agSum[match(names(nclr_per_gene), clr_ds1_agSum[,1]), ]
# Calculate major isoform relative CLR
major_isoform_ratio_clr=clr_ds1_agMax[,2]/clr_ds1_agSum[,2]
# Boxplot
#boxplot(major_isoform_ratio_clr ~ nclr_per_gene, ylim=c(0,1) , main="CLR DS1" , col="lightgrey", ,xlab="N circRNA isoforms per gene", ylab="Major Isoform Relative CLR")
#x=c(1:length(table(nclr_per_gene))); y=1/x; lines(x,y, col="red")
#write.csv(cbind(major_isoform_ratio_clr, nclr_per_gene), "SourceData/Fig3D.MajorIsoform.csv")
#dev.off()

blast=read.csv("../DATA_TABLES/autoBLAST/Output_circ_DS_exonJ_exonJ_whole/stats_circ_DS_exonJ_exonJ_whole_RCM.csv")
colnames(blast)=c("circID.Pair.number", "N_minScore20", "AvgLength_minScore20","N_minScore20_minLength20bp", "AvgCounts_minScore20_minLength20bp", "AvgPercentageAlignment_minScore20", "AvgScore_allMatches", "TotalMatches_Unfiltered" )

rownames(blast)=transpose(strsplit(as.character(blast$circID.Pair.number), split="\t", fixed=TRUE))[[1]]
blast$EnsID=clr_summCTL_ds1$EnsID[match(rownames(blast), rownames(clr_summCTL_ds1))]
na1=which(is.na(blast$EnsID)==TRUE)
blast=blast[-na1, ]

n_per_gene=table(blast$EnsID)
# single=names(n_per_gene)[which(n_per_gene ==1 )]
# blast=blast[-which(blast$EnsID%in%single) , ]

clr_summCTL_ds1_ranked<- clr_summCTL_ds1%>%
  group_by(EnsID) %>%
  mutate(my_ranks = order(order(Mean, decreasing=TRUE)))
clr_summCTL_ds1$Rank=clr_summCTL_ds1_ranked$my_ranks

blast_ranked<- blast %>%
  group_by(EnsID) %>%
  mutate(my_ranks = order(order(N_minScore20, decreasing=TRUE)))
blast$Rank=blast_ranked$my_ranks

clr_summCTL_ds1_use=clr_summCTL_ds1[match(rownames(blast), rownames(clr_summCTL_ds1)) ,]
major=which(clr_summCTL_ds1_use$Rank==1)
major_IsoformData=cbind(clr_summCTL_ds1_use$Coordinate[major], clr_summCTL_ds1_use$EnsID[major], clr_summCTL_ds1_use$Rank[major],  blast$Rank[major])
colnames(major_IsoformData)=c("circId", "EnsID", "clrRank", "blastRank")
major_IsoformData=as.data.frame(major_IsoformData)
major_IsoformData$nIsoforms=n_per_gene[match(major_IsoformData$EnsID, names(n_per_gene))]

major_IsoformData$blastRank=as.numeric(as.character(major_IsoformData$blastRank))
major_IsoformData$nIsoforms=as.numeric(as.character(major_IsoformData$nIsoforms))

t=as.data.frame(table(major_IsoformData$nIsoforms, major_IsoformData$blastRank))
write.csv(table(major_IsoformData$nIsoforms, major_IsoformData$blastRank), "Analysis/DatasetCharacterisationPlots/MajorIsoform_Blast_NminScore20.csv")

# n=10
# plot_blast_results=function(blast, title)
# {
# clr_summCTL_ds1_use=clr_summCTL_ds1[match(rownames(blast), rownames(clr_summCTL_ds1)) ,]
# major=which(clr_summCTL_ds1_use$Rank==1)
# major_IsoformData=cbind(clr_summCTL_ds1_use$Coordinate[major], clr_summCTL_ds1_use$EnsID[major], clr_summCTL_ds1_use$Rank[major],  blast$Rank[major])
# colnames(major_IsoformData)=c("circId", "EnsID", "clrRank", "blastRank")
# major_IsoformData=as.data.frame(major_IsoformData)
# major_IsoformData$nIsoforms=n_per_gene[match(major_IsoformData$EnsID, names(n_per_gene))]
# 
# major_IsoformData$blastRank=as.numeric(as.character(major_IsoformData$blastRank))
# major_IsoformData$nIsoforms=as.numeric(as.character(major_IsoformData$nIsoforms))
# 
# #boxplot(major_IsoformData$blastRank ~ major_IsoformData$nIsoforms)
# t=as.data.frame(table(major_IsoformData$nIsoforms, major_IsoformData$blastRank))
# #write.csv(table(major_IsoformData$nIsoforms, major_IsoformData$blastRank), "/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/Analysis/DatasetCharacterisationPlots/MajorIsoform_Blast_NminScore20.csv")
# 
# library(RColorBrewer)
# coul = brewer.pal(n, "Paired")
# plotdata=t(table(major_IsoformData$nIsoforms, major_IsoformData$blastRank)[1:n,1:n])
# plotdata_percentage=apply(plotdata, 2, function(x){x*100/sum(x,na.rm=T)})
# barplot(plotdata, col=coul, legend.text=TRUE)
# prop.pvals=prop.test(x=plotdata[1,-1], n=apply(plotdata[,-1],2,sum), p=1/c(2:n))$p.value
# barplot(plotdata_percentage, col=coul, sub=prop.pvals, main=title)
# points(100/c(1:n) ~ c(1:n))
# }
# 
# pdf("Analysis/DatasetCharacterisationPlots/autoBlast_Results.pdf", height=11, width=8)
# par(mfrow=c(4,1))
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(N_minScore20, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "N_minScore20")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(AvgLength_minScore20, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "AvgLength_minScore20")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(N_minScore20_minLength20bp, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "N_minScore20_minLength20bp")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(AvgCounts_minScore20_minLength20bp, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "AvgCounts_minScore20_minLength20bp")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(AvgPercentageAlignment_minScore20, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "AvgPercentageAlignment_minScore20")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(AvgScore_allMatches, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "AvgScore_allMatches")
# 
# blast_ranked<- blast %>%
#   group_by(EnsID) %>%
#   mutate(my_ranks = order(order(TotalMatches_Unfiltered, decreasing=TRUE)))
# blast$Rank=blast_ranked$my_ranks
# plot_blast_results(blast, "TotalMatches_Unfiltered")
# 
# dev.off()
