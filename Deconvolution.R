rm(list=ls())
path="PATH_To_PROJECT"
setwd("path")
library(DeconRNASeq)
library(data.table)
library(ggplot2)
library(cowplot)
library(WGCNA)

############################################ LOAD AND FORMAT DATA ############################################

##############Gene count data for DS1 and DS2
load("DATA_TABLES/data_DS1.rda")
load("DATA_TABLES/data_DS2.rda")
geneData_DS1=data_DS1$geneData_DS1$rpkm_filter
circData_DS1=data_DS1$circData_DS1$circCpm_filter

geneData_DS2=data_DS2$geneData_DS2$rpkm_filter
circData_DS2=data_DS2$circData_DS2$circCpm_filter

samples1=read.csv("DATA_TABLES/sampleInfo_DS1.csv")
samples2=read.csv("DATA_TABLES/sampleInfo_DS2.csv")

############## Gencode gene annotation with biotype info 
gencode=read.csv("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GENCODE_human/gencode.v19.annotation.genesOnly.csv")
EnsID=transpose(strsplit(as.character(gencode$gene_id), split=".", fixed=TRUE))
gencode$EnsID=EnsID[[1]]

############## Reference RNA-seq data from the Barres lab (Zhang et al 2016)
#Note: Data from the same cell type are aggregated by "mean"
refB=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/PublishedDatasets/Barres_BrainSeq2.csv")
refB=refB[match(unique(refB[,1]), refB[,1]), ]
rownames(refB)=refB[,1]; refB=refB[,-1]
refB$fetalAstro=apply(refB[,grep("fetal.astrocytes", colnames(refB))], 1,mean)
refB$matureAstro=apply(refB[,grep("mature.astrocytes", colnames(refB))], 1,mean)
refB$Neurons=refB$neurons
refB$Oligo=apply(refB[,grep("Oligodendrocytes", colnames(refB))], 1,mean)
refB$Micro=apply(refB[,grep("Microglia", colnames(refB))], 1,mean)
refB$Endoth=apply(refB[,grep("Endothelial", colnames(refB))], 1,mean)
refB=refB[, 34:39]
#limit to protein coding genes
refB=refB[which(rownames(refB)%in%gencode$gene_name[which(gencode$gene_type%in%"protein_coding")]) ,]

############## Reference RNA-seq data from FANTOM5
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/PublishedDatasets/refF5.pcRNA.rda")
refF5.pcRNA$Symbol=gencode$gene_name[match(rownames(refF5.pcRNA), gencode$EnsID)]
refF5=refF5.pcRNA[-which(is.na(refF5.pcRNA$Symbol)==TRUE) , ]
refF5=refF5[match(unique(refF5$Symbol), refF5$Symbol) , ]
rownames(refF5)=refF5$Symbol; refF5=refF5[, c("Neurons", "Astrocytes")]


#################################################################### DECONVOLVE DS1 and DS2 USING DeconRNAseq
pdf("DeconvolutionPlots_F5vsBarres.pdf", height = 11.7, width = 8.3)
par(mfrow=c(4,4))
par(pch=20, col="grey", family="serif")

###############DS1
data1=geneData_DS1
data1=data1[which(is.na(data1$Symbol)==FALSE) , ]
data1=data1[match(unique(data1$Symbol), data1$Symbol),]
rownames(data1)=data1$Symbol
infocols=c(1,2)
d1 = DeconRNASeq(data1[,-infocols], refB)$out.all
rownames(d1)=colnames(data1)[-infocols]
d1=as.data.frame(d1)

d1F5 = DeconRNASeq(data1[,-infocols], refF5)$out.all
rownames(d1F5)=colnames(data1)[-infocols]
d1F5=as.data.frame(d1F5)

corNP1=cor(d1$Neurons, d1F5$Neurons, method="s")
plot(d1$Neurons, d1F5$Neurons, xlim=c(0,1), ylim=c(0,1), main="DS1 Neuronal Proportion \n refBarres vs refF5", sub=paste0("rho = ", round(corNP1, 3)),
     xlab="refBarres Estimate", ylab="refF5 Estimate")
abline(lm(d1F5$Neurons ~ d1$Neurons), col="red", untf=F)
print(paste0("Correlation between Neurons prop, refBarres vs refF5 in DS1: ", corNP1))
# Conclusion: the correlation between the two deconvolution results is 0.91, so they are expected to give similar results. RefBarres will be used in the following analyses
save(d1, file="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/DATA_TABLES/deconvolution_DS1.rda")
#save(d1F5, file="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/DATA_TABLES/deconvolution_DS1_F5.rda")

###############DS2
data2=geneData_DS2
data2=data2[which(is.na(data2$Symbol)==FALSE) , ]
data2=data2[match(unique(data2$Symbol), data2$Symbol),]
rownames(data2)=data2$Symbol
infocols=c(1,2)
d2 = DeconRNASeq(data2[,-infocols], refB)$out.all
rownames(d2)=colnames(data2)[-infocols]
d2=as.data.frame(d2)

d2F5 = DeconRNASeq(data2[,-infocols], refF5)$out.all
rownames(d2F5)=colnames(data2)[-infocols]
d2F5=as.data.frame(d2F5)

corNP2=cor(d2$Neurons, d2F5$Neurons, method="s")
plot(d2$Neurons, d2F5$Neurons, xlim=c(0,1), ylim=c(0,1), main="DS2 Neuronal Proportion \n refBarres vs refF5", sub=paste0("rho = ", round(corNP2, 3)),
     xlab="refBarres Estimate", ylab="refF5 Estimate")

abline(lm(d2F5$Neurons ~ d2$Neurons), col="red", untf=F)
print(paste0("Correlation between Neurons prop, refBarres vs refF5 in DS2: ", corNP2))

save(d2, file="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/DATA_TABLES/deconvolution_DS2.rda")
# save(d2F5, file="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/DATA_TABLES/deconvolution_DS2_F5.rda")


### SUPP FIG PLOTS
####################################################################  Plot neuronal prop. estimates vs. expression of neuronal-specific genes
pdf("DeconvolutionPlots_MarkerGenes.pdf", height = 11.7, width = 8.3)
par(mfrow=c(2,4))
ng=c("RBFOX1", "MAP2")
for (gene in ng)
{
g=which(data1$Symbol%in%gene)
cor1=cor(d1$Neurons, t(data1[g, -infocols]), method="s")
plot(t(data1[g, -infocols]) ~ d1$Neurons,  main=paste0("DS1: refBarres vs. ", gene), sub=paste0("rho = ", round(cor1,2)), col="indianred",
     xlab = "refBarres Neuronal Proportion", ylab = paste0(gene, " expression (rpkm)"))
cor1=cor(d1F5$Neurons, t(data1[g, -infocols]), method="s")
plot(t(data1[g, -infocols]) ~ d1F5$Neurons, main=paste0("DS1: refF5 vs. ", gene), sub=paste0("rho = ", round(cor1,2)), col="steelblue",
     xlab = "refF5 Neuronal Proportion", ylab = paste0(gene, " expression (rpkm)"))
g=which(data2$Symbol%in%gene)
cor2=cor(d2$Neurons, t(data2[g, -infocols]), method="s")
plot(t(data2[g, -infocols]) ~ d2$Neurons, main=paste0("DS2: refBarres vs. ", gene), sub=paste0("rho = ", round(cor2,2)), col="indianred",
     xlab = "refBarres Neuronal Proportion", ylab = paste0(gene, " expression (rpkm)"))
cor2=cor(d2F5$Neurons, t(data2[g, -infocols]), method="s")
plot(t(data2[g, -infocols]) ~ d2F5$Neurons,  main=paste0("DS2: refF5 vs. ", gene), sub=paste0("rho = ", round(cor2,2)), col="steelblue",
     xlab = "refF5 Neuronal Proportion", ylab = paste0(gene, " expression (rpkm)"))
}
dev.off()


####################################################################  Plot the distribution of neuronal prop. estimates across groups in DS1 and DS2 
pdf("DeconvolutionPlots_ProportionDistributions.pdf", height = 11.7, width = 8.3)

d1$Region=samples1$BroadRegionID[match(rownames(d1), samples1$Sample)]
d1$Pheno=samples1$ASD.CTL[match(rownames(d1), samples1$Sample)]
d1$Region_Pheno=paste(d1$Region, d1$Pheno, sep="_")

d1F5$Region=samples1$BroadRegionID[match(rownames(d1F5), samples1$Sample)]
d1F5$Pheno=samples1$ASD.CTL[match(rownames(d1F5), samples1$Sample)]
d1F5$Region_Pheno=paste(d1F5$Region, d1F5$Pheno, sep="_")

d2$Region=samples2$BroadRegionID[match(rownames(d2), samples2$Sample)]
d2$Pheno=samples2$ASD.CTL[match(rownames(d2), samples2$Sample)]
d2$Region_Pheno=paste(d2$Region, d2$Pheno, sep="-")

d2F5$Region=samples2$BroadRegionID[match(rownames(d2F5), samples2$Sample)]
d2F5$Pheno=samples2$ASD.CTL[match(rownames(d2F5), samples2$Sample)]
d2F5$Region_Pheno=paste(d2F5$Region, d2F5$Pheno, sep="_")

p1 <- ggplot(as.data.frame(d1) , aes(y=Neurons, x=Region_Pheno, fill=Region_Pheno)) +
  geom_point() +  geom_violin(alpha=0.2) + geom_boxplot(width=0.1, fill="white" ,alpha=0.2) +
  ggtitle("DS1, refBarres") + ylab("Neuronal Proportion (refBarres)") + theme_bw() + 
  theme(legend.position = "null") + xlab("Region & Phenotype")

p2 <- ggplot(as.data.frame(d2) , aes(y=Neurons, x=Region_Pheno, fill=Region_Pheno)) +
  geom_point() +  geom_violin(alpha=0.2) + geom_boxplot(width=0.1, fill="white" ,alpha=0.2) +
  ggtitle("DS2, refBarres") + ylab("Neuronal Proportion (refBarres)") + theme_bw() + 
  theme(legend.position = "null") + xlab("Region & Phenotype")

p3 <- ggplot(as.data.frame(d1F5) , aes(y=Neurons, x=Region_Pheno, fill=Region_Pheno)) +
  geom_point() +  geom_violin(alpha=0.2) + geom_boxplot(width=0.1, fill="white" ,alpha=0.2) +
  ggtitle("DS1, refF5") + ylab("Neuronal Proportion (refF5)") + theme_bw() + 
  theme(legend.position = "null") + xlab("Region & Phenotype")

p4 <- ggplot(as.data.frame(d2F5) , aes(y=Neurons, x=Region_Pheno, fill=Region_Pheno)) +
  geom_point() +  geom_violin(alpha=0.2) + geom_boxplot(width=0.1, fill="white" ,alpha=0.2) +
  ggtitle("DS2, refF5") + ylab("Neuronal Proportion (refF5)") + theme_bw() + 
  theme(legend.position = "null") + xlab("Region & Phenotype")

p5 <- ggplot(as.data.frame(d2F5) , aes(y=Neurons, x=Region_Pheno, fill=Region_Pheno)) +
  geom_blank()

plot_grid(p1, p2, p3, p4, p5, p5, nrow = 3, ncol = 2)
dev.off()


#################################################################### Neuronal prop. estimates vs gene exp PCA
pdf("DeconvolutionPlots_PCAbyGenes.pdf", height = 11.7, width = 8.3)
xlab="Neuronal Proportion (refBarres)"
ylab="PC1 of Gene Expression"
par(mfrow=c(2,2))
samples1$colors <- factor(samples1$RegionID, levels=(c("ba9", "ba41-42-22", "vermis")) ,labels=c("#66c2a5",  "#8da0cb" ,"#fc8d62"))
samples1$pch=17;samples1$pch[which(samples1$ASD.CTL %in%"ASD")]=16
samples2$colors <- factor(samples2$RegionID, levels=(c("ba9", "ba41-42-22", "vermis")) ,labels=c("#66c2a5",  "#8da0cb" ,"#fc8d62"))
samples2$pch=17;samples2$pch[which(samples2$ASD.CTL %in%"ASD")]=16

###All regions
data1=data_DS1$geneData_DS1$rpkm_filter; data2=data_DS2$geneData_DS2$rpkm_filter
pca1=princomp(cor(data1[, -c(1,2)], method="s"))$loadings
pca2=princomp(cor(data2[, -c(1,2)], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1], col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: All Regions")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2], col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: All Regions")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")

###CTX
pca1=princomp(cor(data1[, which(colnames(data1)%in%samples1$Sample[which(samples1$BroadRegionID=="CTX")])], method="s"))$loadings
pca2=princomp(cor(data2[, which(colnames(data2)%in%samples2$Sample[which(samples2$BroadRegionID=="CTX")])], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1], col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: CTX")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2], col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: CTX")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")

###CB
pca1=princomp(cor(data1[, which(colnames(data1)%in%samples1$Sample[which(samples1$BroadRegionID=="CB")])], method="s"))$loadings
pca2=princomp(cor(data2[, which(colnames(data2)%in%samples2$Sample[which(samples2$BroadRegionID=="CB")])], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1],col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: CB")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2],col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: CB")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")
dev.off()


#################################################################### Neuronal prop. estimates vs circRNA PCA
pdf("DeconvolutionPlots_PCAbyCirc.pdf", height = 11.7, width = 8.3)
par(mfrow=c(4,2))
xlab="Neuronal Proportion (refBarres)"
ylab="PC1 of circRNA Expression"
rm(data1); rm(data2)
data1=data_DS1$circData_DS1$circCpm
data2=data_DS2$circData_DS2$circCpm

###All regions
pca1=princomp(cor(data1[, -c(1:15)], method="s"))$loadings
pca2=princomp(cor(data2[, -c(1:15)], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1], col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: All Regions")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2], col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: All Regions")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")


###CTX
pca1=princomp(cor(data1[, which(colnames(data1)%in%samples1$Sample[which(samples1$BroadRegionID=="CTX")])], method="s"))$loadings
pca2=princomp(cor(data2[, which(colnames(data2)%in%samples2$Sample[which(samples2$BroadRegionID=="CTX")])], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1], col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: CTX")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2], col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: CTX")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")


###CB
pca1=princomp(cor(data1[, which(colnames(data1)%in%samples1$Sample[which(samples1$BroadRegionID=="CB")])], method="s"))$loadings
pca2=princomp(cor(data2[, which(colnames(data2)%in%samples2$Sample[which(samples2$BroadRegionID=="CB")])], method="s"))$loadings
m1=match(rownames(pca1), rownames(d1)); colors1=as.character(samples1$colors[match(rownames(pca1), samples1$Sample)]); pch1=samples1$pch[match(rownames(pca1), samples1$Sample)]
m2=match(rownames(pca2), rownames(d2)); colors2=as.character(samples2$colors[match(rownames(pca2), samples2$Sample)]); pch2=samples2$pch[match(rownames(pca2), samples2$Sample)]

plot(pca1[,1] ~ d1$Neurons[m1],col=colors1,pch=pch1, sub=paste0("rho = ",round(cor(pca1[,1] , d1$Neurons[m1], method="s"),2)), ylab=ylab, xlab=xlab, main="DS1: CB")
abline(lm(pca1[,1] ~ d1$Neurons[m1]), col="grey")
plot(pca2[,1] ~ d2$Neurons[m2],col=colors2,pch=pch2, sub=paste0("rho = ",round(cor(pca2[,1] , d2$Neurons[m2], method="s"),2)), ylab=ylab, xlab=xlab, main="DS2: CB")
abline(lm(pca2[,1] ~ d2$Neurons[m2]), col="grey")

dev.off()

