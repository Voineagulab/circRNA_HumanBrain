setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
rm(list=ls())
library(WGCNA)
library(data.table)
library(splines)
library(ggplot2)
###################################################################################################
makeSTARGeneCountsTable <- function(path, seqType){
  files <- read.table(file = paste0(path, "/", "STARSummary.txt"), check.names=FALSE, as.is = TRUE, sep = "\t", fill = TRUE, header=TRUE)
  ### generate gene counts table
  for (j in c(2: ncol(files))){
    sample <- as.character(colnames(files)[j])
    print(sample)
    data <- read.table(paste0(path, "/", sample, "/", sample, "ReadsPerGene.out.tab"), sep="\t")
    data <- data[-c(1:4),]
    colnames(data) <- c("EnsID", "Unstranded", "First", "Second")
    if (j==2){
      geneCounts=data[, c("EnsID",seqType)]
    }else{
      m <- match(geneCounts$EnsID, data$EnsID)
      geneCounts <- cbind(geneCounts, data[m, seqType])
    }
    colnames(geneCounts)[ncol(geneCounts)] <- sample
  }
  rownames(geneCounts) <- geneCounts$EnsID
  
  return(geneCounts)
}
###################################################################################################
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
############################################ LOAD DATA ############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
#Extract filtered data for genes, circRNAs and splice junctions,  DS1
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
#Extract filtered data for genes, circRNAs and splice junctions,  DS2
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
# Define infocols
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)
############################################ Neur vs Astro, our data ##############################
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/Analysis/Neurons_Astro/circCountAstroFilter.rda")
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/Analysis/Neurons_Astro/circCountNeuronFilter.rda")
neur_annot=read.csv("Analysis/Neurons_Astro/Sample_Annotation.csv")
#use astrocyte cpm data. Select one of the astrocyte samples, so that N is the same as for neurons
astro=circCountAstroFilter.df[rowSums(circCountAstroFilter.df[, grep("cpm", colnames(circCountAstroFilter.df))] >=0.1) >=1 , ]
colSums(astro[, grep("cpm", colnames(astro))] >=0.1) 
#rbVect1_S1_cpm rbVect2_S2_cpm rbVect3_S3_cpm 
#89             93            107 
# Will report the middle value: N=93
astro=astro[which(astro$rbVect2_S2_cpm >=0.1) , ]
# use control neurons cpm data
neur=circCountNeuronFilter.df[which(circCountNeuronFilter.df[,grep("GOK3380A1", colnames(circCountNeuronFilter.df))[2]]>=0.1), ]
astro$Id=paste(paste0("chr", astro$Chr), astro$Start, astro$End, sep="_")
neur$Id=paste(paste0("chr",neur$Chr), neur$Start, neur$End, sep="_")

dim(neur); dim(astro); length(intersect(neur$Id, astro$Id)) 
#[1] 3601   neuronal circRNAs
#[1] 93  astrocyte circRNAs
#[1] 74 common circRNAs

#rownames(neur)=paste0("chr", gsub("-", "_", neur$Coordinate, fixed=TRUE))
#Overlap with DS1
circ_ds1=data_DS1$circData_DS1$circCounts_filter
length(intersect(neur$Id, circ_ds1$Id))
#2913
length(intersect(astro$Id, circ_ds1$Id))
#90

#> 90/93
# [1] 0.9677419
# > 2913/3601
# [1] 0.808942

neur_bed=neur[, c(2,3,4,9,8,7)]
neur_bed[,5]=0
neur_bed[,6]="."
neur_bed[,1]=paste0("chr", neur_bed[,1])

astro_bed=astro[, c(2,3,4,13,8,7)]
astro_bed[,5]=0
astro_bed[,6]="."
astro_bed[,1]=paste0("chr", astro_bed[,1])

write.table(neur_bed, "Analysis/neur.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(astro_bed, "Analysis/astro.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

neur_mm10=read.table("Analysis/neur_mm10_lift.bed", sep="\t")
astro_mm10=read.table("Analysis/astro_mm10_lift.bed", sep="\t")
mm10_nuclear=read.csv("Analysis/Mouse_vs_Human/NIHMS763145-supplement.csv")

astro_mm10$Id=paste0(astro_mm10[,1], ":", astro_mm10[,2] , "-", astro_mm10[,3]+1)
neur_mm10$Id=paste0(neur_mm10[,1], ":", neur_mm10[,2], "-", neur_mm10[,3]+1)

astro_mm10=cbind(astro_mm10, astro[match(astro_mm10[,4], astro$Id) ,c("Gene", "rbVect2_S2_cpm")])

length(intersect(mm10_nuclear$Location, astro_mm10$Id))
#[1] 10
length(intersect(mm10_nuclear$Location, neur_mm10$Id))
#[1] 100
10/nrow(astro_mm10)
#[1] 0.1075269
100/nrow(neur_mm10)
#[1]0.02884338
colSums(nuclear[, grep("cpm", colnames(mm10_nuclear))] >0)
# cpm.N cpm.A cpm.O cpm.C 
# 2590  1512  1733  2034 
colSums(nuclear[, grep("cpm", colnames(mm10_nuclear))] > 0.1)
# cpm.N cpm.A cpm.O cpm.C 
# 158    65    37  2034 
############################################# Pasca lab data from brain organoids ################
############## #LOAD DATA ###############
#setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/")
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/Sloan_Pasca_DATA/DATA_TABLES/circCoordinateAndCpm.rda")
#load("/Volumes/Data1/PROJECTS/circRNAs/DHG_DATA/PublishedDatasets/Steven_Solan_Neuron_data/DATA_TABLES/circCoordinateAndCount.rda")

infocols=c(1:9)
plab_cpm=circCoordinateAndCpm.df
##Filter
plab_cpm=plab_cpm[rowSums(plab_cpm[, -infocols] > 0) >= 2 , ]
plab_cpm_exp=plab_cpm[,-infocols];rownames(plab_cpm_exp)=paste0("chr", plab_cpm$Coordinate)
dim(plab_cpm)
#[1] 5023   31
##Overlap with DS1
plab_cpm_id=paste(paste0("chr", plab_cpm$Chr), plab_cpm$Start, plab_cpm$End, sep="_")

overlap=which(plab_cpm_id%in% unique(c(data_DS1$circData_DS1$circCpm_filter$Id , data_DS2$circData_DS2$circCpm_filter$Id)))
length(overlap)
#[1] 2820
length(overlap)/nrow(plab_cpm)
#0.5614175

# #Checking whether the libsize has an effect on NcircExp after CPM norm
# plab_libsizes=plab_libsizes[match(colnames(plab_cpm)[-infocols], names(plab_libsizes))]
# plab_libsizes=as.numeric(as.character(plab_libsizes))
# names(plab_libsizes)=colnames(plab_counts)[-infocols]
# plab_nExp=colSums(plab_cpm[,-infocols] >= 0.1)
# plot(plab_libsizes/1000000,plab_nExp)
# # looks ok.

##########################################  Number of expressed circRNAs vs maturation timepoint ############## 
############## get timepoint
tp=transpose(strsplit(colnames(plab_cpm_exp), "_"))[[5]]
tp=gsub("b", "", tp); tp=as.numeric(tp)
tp_bin=cut(tp, seq(0,600,200))
############## select neuronal and astrocyte (i.e. glial) samples 
plab_neurons=grep("Thy", colnames(plab_cpm_exp))
plab_glia=grep("Hepa", colnames(plab_cpm_exp))
############## identify circRNAs expressed in neurons and astrocytes (i.e. glia) 
plab_NeuronsCirc=which(rowSums(plab_cpm_exp[,intersect(which(tp>=150) ,grep("Thy", colnames(plab_cpm_exp))) ] >=0.1) >=2)
plab_GliaCirc=which(rowSums(plab_cpm_exp[,intersect(which(tp>=150) ,grep("Hepa", colnames(plab_cpm_exp))) ] >=0.1) >=2)
length(plab_NeuronsCirc)
#[1] 420
length(plab_GliaCirc)
#[1] 318
############## Number of circRNAs expressed at a minimum of 0.1 CPM in each sample
nCircExp=colSums(plab_cpm_exp[,] >= 0.1)  
############## Generate a dataframe for plotting in ggplot2
d=as.data.frame(cbind(nCircExp, tp, "astrocytes"))
colnames(d)=c("nCircExp", "TimePoint", "CellType")
d$CellType=as.character(d$CellType); d$CellType[plab_neurons]="neurons"; 
d$nCircExp=as.numeric(as.character(d$nCircExp)) 
d$TimePoint=as.numeric(as.character(d$TimePoint))
##############Plot                                                                          
pdf("Analysis/Neurons_Astro/Fig6.nCircExp_vs_Timepoint.pdf", height=4, width=6)
p1=ggplot(d, aes(x=TimePoint, y=nCircExp, color=CellType)) +
  geom_point(size=5, shape=17) + 
  geom_smooth(method=lm, linetype=2) +
  scale_colour_manual(values=c("plum3","steelblue3" ))
print(p1)
dev.off()
##############Linear model for N circ vs timepoint in mature cells
mature=which(d$TimePoint >= 150)
summary(lm(d$nCircExp[mature] ~ d$TimePoint[mature] + d$CellType[mature]))
write.csv(d, "SourceData/Fig7A.csv")
# Call:
#   lm(formula = d$nCircExp[mature] ~ d$TimePoint[mature] + d$CellType[mature])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -122.47  -51.37  -32.84   40.40  205.84 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                606.199     88.341   6.862 7.37e-05 ***
#   d$TimePoint[mature]         -1.067      0.266  -4.010  0.00306 ** 
#   d$CellType[mature]neurons  226.357     68.256   3.316  0.00899 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 109.9 on 9 degrees of freedom
# Multiple R-squared:  0.7831,	Adjusted R-squared:  0.7349 
# F-statistic: 16.24 on 2 and 9 DF,  p-value: 0.001031
############# Linear model for individual circRNAs in neurons and glia
tp_n=tp[plab_neurons]
tp_g=tp[plab_glia]

modelN=lm(t(plab_cpm_exp[, plab_neurons]) ~ tp_n)
pvalsN=adjPvals.from.lm(coefficients(summary(modelN)))[,2]
corPN=cor(t(plab_cpm_exp[, plab_neurons]), tp_n, method="p")
corSN=cor(t(plab_cpm_exp[, plab_neurons]), tp_n, method="s")


modelG=lm(t(plab_cpm_exp[, plab_glia]) ~ tp_g)
pvalsG=adjPvals.from.lm(coefficients(summary(modelG)))
corPG=cor(t(plab_cpm_exp[, plab_glia]), tp_g, method="p")
corSG=cor(t(plab_cpm_exp[, plab_glia]), tp_g, method="s")

DE=data.frame(pvalsN, corPN, corSN, pvalsG, corPG, corSG)
rownames(DE)=rownames(plab_cpm_exp)
############# Gene Exp
# Counts
path="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/Sloan_Pasca_DATA/RESULTS/STAR_Output/"
seqType="Unstranded"
plab_GeneCounts=makeSTARGeneCountsTable(path, seqType)
# CPM
summary=read.table(header=TRUE,"/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/Sloan_Pasca_DATA/STARSummary.txt")
plab_libsizes=summary[grep("Uniquely mapped reads number |", summary$summary, fixed=TRUE), -1]
plab_libsizes=unlist(plab_libsizes)
snames=gsub(".", "-", names(plab_libsizes) ,fixed=TRUE)
plab_libsizes=as.numeric(as.character(plab_libsizes))
names(plab_libsizes)=snames
m=match(colnames(plab_GeneCounts), names(plab_libsizes))
plab_GeneCpm=plab_GeneCounts
for (j in c(2:ncol(plab_GeneCounts))) plab_GeneCpm[,j]=plab_GeneCounts[,j]/(plab_libsizes[m[j]]/10^6)
#RPKM
load("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GENCODE_human/exonicLength.rda")
rownames(exonicLength)=transpose(strsplit(rownames(exonicLength), ".", fixed=TRUE))[[1]]
length=transpose(exonicLength)[[1]]; names(length)=rownames(exonicLength)
length=length[match(rownames(plab_GeneCpm), names(length))]/1000
plab_RPKM=plab_GeneCpm
for (j in c(2:ncol(plab_RPKM))) plab_RPKM[,j]=plab_GeneCpm[,j]/length
plab_RPKM=plab_RPKM[-which(is.na(length)==TRUE),]
#Filter
plab_RPKM=plab_RPKM[ rowSums(plab_RPKM[,-1] >= 1) >= ncol(plab_RPKM)/2 , ]
dim(plab_RPKM)
overlap_G=intersect(data_DS1$geneData_DS1$rpkm_filter$EnsID , plab_RPKM$EnsID)
length(overlap_G)
length(overlap_G)/nrow(plab_RPKM)

#######Save
data_plab=list(circ=plab_cpm, genes=plab_RPKM)
save(data_plab, file="DATA_TABLES/data_plab.rda")
############## Generate fasta files for MEME-ChIP
##############Astrocyte-specific circRNAs
astro_circ=rownames(plab_cpm_exp)[setdiff(plab_GliaCirc, plab_NeuronsCirc)]
#astro_circ=rownames(plab_cpm_exp)[plab_GliaCirc]
#astro_circ=rownames(plab_cpm_exp)[intersect(plab_GliaCirc, which(plab_cpm_exp[, grep("Thy1_Day_450", colnames(plab_cpm_exp))]==0))]
astro_circ_bed="/Analysis/DE/astro_circ.bed"
astro_circ_fasta="/Analysis/DE/astro_circ.fasta"
generate.fasta_for_MEMEChIP(astro_circ, w=100, astro_circ_bed, astro_circ_fasta)

##############Neuron-specific circRNAs
neur_circ=rownames(plab_cpm_exp)[setdiff(plab_NeuronsCirc, plab_GliaCirc)]
#neur_circ=rownames(plab_cpm_exp)[intersect(plab_NeuronsCirc, which(plab_cpm_exp[, grep("Hepa_Day_450", colnames(plab_cpm_exp))]==0))]
neur_circ_bed="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/Analysis/DE/neur_circ.bed"
neur_circ_fasta="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/Analysis/DE/neur_circ.fasta"
generate.fasta_for_MEMEChIP(neur_circ, w=100, neur_circ_bed, neur_circ_fasta)

##############Control (i.e. background) sequences
exons=read.table("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GTF/GTF_to_BED/exonsUnique.bed", sep="\t", header=FALSE)
gencode=read.csv("/Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/ANNOTATION/GENCODE_human/gencode.v19.annotation.genesOnly.csv")
EnsID=transpose(strsplit(as.character(gencode$gene_id), split=".", fixed=TRUE))
gencode$EnsID=EnsID[[1]]
elim=unique(c(data_DS1$circData_DS1$circCpm_filter$Symbol, data_DS2$circData_DS2$circCpm_filter$Symbol, gencode$EnsID[match(plab_cpm$Gene, gencode$gene_name)]))
exons=exons[-which(exons$V5%in%elim), ]
exons=exons[which(exons$V5%in%data_DS1$geneData_DS1$rpkm_filter$EnsID), ]
s=sample(c(1:nrow(exons)), size=100000)
control_seq=gsub("exon_", "", exons$V4[s], fixed=TRUE)
control_bed="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/Analysis/DE/meme_control.bed"
control_fasta="/Volumes/Data0/PROJECTS/circRNAs/DHG_DATA/CompleteDataset/Analysis/DE/meme_control.fasta"    
generate.fasta_for_MEMEChIP(control_seq, w=100, control_bed, control_fasta)

