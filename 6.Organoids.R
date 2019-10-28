
rm(list=ls())
library(WGCNA)
library(data.table)
library(splines)
library(ggplot2)

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

############################################# Pasca lab data from brain organoids ################
############## #LOAD DATA ###############
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/Sloan_Pasca_DATA/DATA_TABLES/circCoordinateAndCpm.rda")
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
