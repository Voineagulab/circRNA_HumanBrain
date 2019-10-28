rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
se=read.table("../DATA_TABLES/rMATS_OUTPUT/MATS_output/SE.MATS.JunctionCountOnly.txt", sep="\t", header=TRUE)
se$chr=paste0("chr", se$chr)
#ds=which ((se$FDR < 0.05) & (abs(se$IncLevelDifference) > 0.1))
ds=which ((se$FDR < 0.05) & (abs(se$IncLevelDifference) > 0))
se=se[ds,]

se$startCoord=paste(se$chr, se$exonStart_0base +1, se$strand, sep="_")
se$endCoord=paste(se$chr, se$exonEnd, se$strand, sep="_")

de=read.csv("Analysis/DE/SuppTable.region_circ.csv")

de$startCoord=paste(de$Chr, de$Start, de$Strand, sep="_")
de$endCoord=paste(de$Chr, de$End, de$Strand, sep="_")

de_ds=which((de$startCoord%in%se$startCoord)|(de$endCoord%in%se$endCoord))
