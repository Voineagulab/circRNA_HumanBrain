rm(list=ls())
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/DATA_TABLES/")
source('/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/SCRIPTS/Final/Supplementary_Code/Functions.R')
#################################################################### LOAD DATA
###### Read in data from DS1 and DS2 data folders
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch2/RESULTS_DS1/DATA_TABLES_DS1/circData_DS1.rda")
circData_DS1=circData
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch1/RESULTS_DS2/DATA_TABLES_DS2/circData_DS2.rda")
circData_DS2=circData
rm(circData)
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch2/RESULTS_DS1/DATA_TABLES_DS1/geneData_DS1.rda")
geneData_DS1=geneData
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch1/RESULTS_DS2/DATA_TABLES_DS2/geneData_DS2.rda")
geneData_DS2=geneData
rm(geneData)
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch2/RESULTS_DS1/DATA_TABLES_DS1/sjData_DS1.rda")
sjData_DS1=sjData
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch1/RESULTS_DS2/DATA_TABLES_DS2/sjData_DS2.rda")
sjData_DS2=sjData
rm(sjData)
###### Read in Sample Info
sampleInfo_DS1=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch2/RESULTS_DS1/sampleNames_Info_complete_DS1.csv")
sampleInfo_DS2=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch1/RESULTS_DS2/sampleNames_Info_complete_DS2.csv")

#################################################################### Format Circ RNA data
####Add "chr" to the chromosome column and add a "Coordinate" column that includes strand information. There will be 15 infocols after this.
for (j in c(1:length(circData_DS1)))
{
circData_DS1[[j]]$Id <- sub("^", "chr", circData_DS1[[j]]$Id)
circData_DS1[[j]]$Chr <- sub("^", "chr", circData_DS1[[j]]$Chr)
Coordinate=paste(circData_DS1[[j]]$Id, circData_DS1[[j]]$Strand, sep="_")
circData_DS1[[j]]=as.data.frame(cbind(Coordinate, circData_DS1[[j]]))
circData_DS1[[j]][, c(1:5)]=circData_DS1[[j]][, c("Coordinate","Chr","Start","End", "Id")]
colnames(circData_DS1[[j]])[c(1:5)]=c("Coordinate","Chr","Start","End", "Id")

circData_DS2[[j]]$Id <- sub("^", "chr", circData_DS2[[j]]$Id)
circData_DS2[[j]]$Chr <- sub("^", "chr", circData_DS2[[j]]$Chr)
Coordinate=paste(circData_DS2[[j]]$Id, circData_DS2[[j]]$Strand, sep="_")
circData_DS2[[j]]=as.data.frame(cbind(Coordinate, circData_DS2[[j]]))

rownames(circData_DS1[[j]])=circData_DS1[[j]]$Coordinate
rownames(circData_DS2[[j]])=circData_DS2[[j]]$Coordinate
}
#################################################################### Format SJ data
####Add "chr" to the chromosome column and add a "Coordinate" column that includes strand information. There will be 15 infocols after this.
for (j in c(1:length(sjData_DS1)))
{ 
  sjData_DS1[[j]]$Id <- sub("^", "chr", sjData_DS1[[j]]$Id)
  sjData_DS1[[j]]$Chr <- sub("^", "chr", sjData_DS1[[j]]$Chr)
  Coordinate=paste(sjData_DS1[[j]]$Id, sjData_DS1[[j]]$Strand, sep="_")
  sjData_DS1[[j]]=as.data.frame(cbind(Coordinate, sjData_DS1[[j]]))
  sjData_DS1[[j]][, c(1:5)]=sjData_DS1[[j]][, c("Coordinate","Chr","Start","End", "Id")]
  colnames(sjData_DS1[[j]])[c(1:5)]=c("Coordinate","Chr","Start","End", "Id")
  
  sjData_DS2[[j]]$Id <- sub("^", "chr", sjData_DS2[[j]]$Id)
  sjData_DS2[[j]]$Chr <- sub("^", "chr", sjData_DS2[[j]]$Chr)
  Coordinate=paste(sjData_DS2[[j]]$Id, sjData_DS2[[j]]$Strand, sep="_")
  sjData_DS2[[j]]=as.data.frame(cbind(Coordinate, sjData_DS2[[j]]))
  
  rownames(sjData_DS1[[j]])=sjData_DS1[[j]]$Coordinate
  rownames(sjData_DS2[[j]])=sjData_DS2[[j]]$Coordinate
  
}
#################################################################### Format Sample Info files and libsizes 
#### Remove duplicated columns from sample info, and the samples that have been eliminated at QC level
sampleInfo_DS1=sampleInfo_DS1[, -which(colnames(sampleInfo_DS1)%in%c("BrainID.1","Region","ASD.CTL.1"))]
sampleInfo_DS1=sampleInfo_DS1[sampleInfo_DS1$ElimQC=="No" , ];sampleInfo_DS2=sampleInfo_DS2[sampleInfo_DS2$ElimQC=="No" , ]
#### Add cellular composition info to Sample Info
load("deconvolution_DS1.rda")
load("deconvolution_DS1_F5.rda")
colnames(d1F5)=paste("F5", colnames(d1F5), sep="_")
rownames(sampleInfo_DS1)=sampleInfo_DS1$Sample
sampleInfo_DS1=merge(sampleInfo_DS1, d1, by="row.names")
rownames(sampleInfo_DS1)=sampleInfo_DS1$Sample; sampleInfo_DS1=sampleInfo_DS1[,-1]
sampleInfo_DS1=merge(sampleInfo_DS1, d1F5, by="row.names")
rownames(sampleInfo_DS1)=sampleInfo_DS1$Sample; sampleInfo_DS1=sampleInfo_DS1[,-1]

load("deconvolution_DS2.rda")
load("deconvolution_DS2_F5.rda")
colnames(d2F5)=paste("F5", colnames(d2F5), sep="_")

rownames(sampleInfo_DS2)=sampleInfo_DS2$Sample

sampleInfo_DS2=merge(sampleInfo_DS2, d2, by="row.names")
rownames(sampleInfo_DS2)=sampleInfo_DS2$Sample; sampleInfo_DS2=sampleInfo_DS2[,-1]

sampleInfo_DS2=merge(sampleInfo_DS2, d2F5, by="row.names")
rownames(sampleInfo_DS2)=sampleInfo_DS2$Sample; sampleInfo_DS2=sampleInfo_DS2[,-1]
#### Add "Broad region" info to the sample info file
sampleInfo_DS1$BroadRegionID=NA
sampleInfo_DS1$BroadRegionID=as.character(sampleInfo_DS1$BroadRegionID)
sampleInfo_DS1$BroadRegionID[grep("ba", sampleInfo_DS1$RegionID)]="CTX"
sampleInfo_DS1$BroadRegionID[-grep("ba", sampleInfo_DS1$RegionID)]="CB"

sampleInfo_DS2$BroadRegionID=NA
sampleInfo_DS2$BroadRegionID=as.character(sampleInfo_DS2$BroadRegionID)
sampleInfo_DS2$BroadRegionID[grep("ba", sampleInfo_DS2$RegionID)]="CTX"
sampleInfo_DS2$BroadRegionID[-grep("ba", sampleInfo_DS2$RegionID)]="CB"

##### Format libsizes
summary1=read.table("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch2/RESULTS_DS1/STAR_GeneCounts_DS1/STARSummary.txt", header=TRUE)
libsizes1=summaryToLibsizes(summary1)

summary2=read.table("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/Batch1/RESULTS_DS2/STAR_GeneCounts_DS2/summaryTable.txt", header=TRUE)
libsizes2=summaryToLibsizes(summary2)

#################################################################### Select the filtering threshold for CircRNA data
ns1=c(1:10); ns2=c(1:10)
resultFilterCounts=matrix(NA, ncol=4, nrow=length(ns1)); colnames(resultFilterCounts)=c("DS1", "DS2", "Common", "Common/DS2"); rownames(resultFilterCounts)=paste("2counts_",ns1, "DS1", ns2, "DS2", sep="_")
resultFilterCpm=matrix(NA, ncol=4, nrow=length(ns1)); colnames(resultFilterCpm)=c("DS1", "DS2", "Common" ,"Common/DS2"); rownames(resultFilterCpm)=paste("0.1cpm_", ns1, "DS1", ns2, "DS2", sep="_")
infocols=c(1:15)

for (s in c(1:length(ns1))) 
{ 
th1=ns1[s]
th2=ns2[s]
ds1=circData_DS1$circCounts[rowSums(circData_DS1$circCounts[, -infocols] >= 2) >= th1 , ]; resultFilterCounts[s,1]=nrow(ds1)
ds2=circData_DS2$circCounts[rowSums(circData_DS2$circCounts[, -infocols] >= 2) >= th2 , ]; resultFilterCounts[s,2]=nrow(ds2)
resultFilterCounts[s,3]=length(intersect(ds1$Coordinate, ds2$Coordinate))
resultFilterCounts[s,4]=round(length(intersect(ds1$Coordinate, ds2$Coordinate))/nrow(ds2), 3)

ds1=circData_DS1$circCounts[rowSums(circData_DS1$circCpm[, -infocols] >= 0.1) >= th1 , ]; resultFilterCpm[s,1]=nrow(ds1)
ds2=circData_DS2$circCounts[rowSums(circData_DS2$circCpm[, -infocols] >= 0.1) >= th2 , ]; resultFilterCpm[s,2]=nrow(ds2)
resultFilterCpm[s,3]=length(intersect(ds1$Coordinate, ds2$Coordinate))
resultFilterCpm[s,4]=round(length(intersect(ds1$Coordinate, ds2$Coordinate))/nrow(ds2), 3)
}
result=rbind(resultFilterCounts, resultFilterCpm)
write.csv(result, "../Analysis/DatasetCharacterisationPlots/SuppFig.filteringComparisons_all_ds1_ds2.csv")
##### Based on the above, 2 counts in 5 samples seems like a reasonable permissive filter
circData_DS1$circCounts_filter=circData_DS1$circCounts[rowSums(circData_DS1$circCounts[, -infocols] >= 2) >= 5 , ]
circData_DS1$circCpm_filter=circData_DS1$circCpm[which(rownames(circData_DS1$circCpm)%in%rownames(circData_DS1$circCounts_filter)) ,]
circData_DS2$circCounts_filter=circData_DS2$circCounts[rowSums(circData_DS2$circCounts[, -infocols] >= 2) >= 5 , ]
circData_DS2$circCpm_filter=circData_DS2$circCpm[which(rownames(circData_DS2$circCpm)%in%rownames(circData_DS2$circCounts_filter)) ,]

#################################################################### Filter the sjData based on the filtered circRNAs
for (j in c(1:4))
{
k=j+4
sjData_DS1[[k]]=sjData_DS1[[j]][which(rownames(sjData_DS1[[j]])%in%rownames(circData_DS1$circCounts_filter)),]
names(sjData_DS1)[[k]]=paste(names(sjData_DS1)[[j]], "filter", sep="_")
sjData_DS2[[k]]=sjData_DS2[[j]][which(rownames(sjData_DS2[[j]])%in%rownames(circData_DS2$circCounts_filter)),]
names(sjData_DS2)[[k]]=paste(names(sjData_DS2)[[j]], "filter", sep="_")
}
#################################################################### Put all data with columns in the same order as the sampleinfo rows

circData_DS1=reorder.f(circData_DS1, sampleInfo_DS1, infocols)
circData_DS2=reorder.f(circData_DS2, sampleInfo_DS2, infocols)
sjData_DS1=reorder.f(sjData_DS1, sampleInfo_DS1, infocols)
sjData_DS2=reorder.f(sjData_DS2, sampleInfo_DS2, infocols)
geneData_DS1=reorder.f(geneData_DS1, sampleInfo_DS1, c(1,2))
geneData_DS2=reorder.f(geneData_DS2, sampleInfo_DS2, c(1,2))

libsizes1=libsizes1[match(rownames(sampleInfo_DS1), names(libsizes1))]
libsizes2=libsizes2[match(rownames(sampleInfo_DS2), names(libsizes2))]

#################################################################### Generate Lists
data_DS1=list(circData_DS1=circData_DS1, sjData_DS1=sjData_DS1,geneData_DS1=geneData_DS1,sampleInfo_DS1=sampleInfo_DS1,libsizes_DS1=libsizes1)
data_DS2=list(circData_DS2=circData_DS2, sjData_DS2=sjData_DS2,geneData_DS2=geneData_DS2,sampleInfo_DS2=sampleInfo_DS2,libsizes_DS2=libsizes2)
#################################################################### UpdateGene Symbols
t_s=read.table("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/ensemblToGeneName_02.2019.txt", sep="\t", header=FALSE)
t_g=read.table("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/ensGene.02.2019.txt", sep="\t", header=FALSE)
ens=as.data.frame(cbind(t_s, t_g$V13[match(t_s$V1, t_g$V2)]))
colnames(ens)=c("TxID", "Symbol", "EnsID")
save(ens,file="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/ensGeneTxSymbol.02.2019.rda")

getSymbols=function(x, ref) s=ref$Symbol[match(x, ref$EnsID)]
for (j in c("circData_DS1"  , "sjData_DS1"   ,  "geneData_DS1"))
  for (k in c(1: length(data_DS1[[j]])))
  { print(paste(j, "____", k))
    print("Before:")
    print(length(which(data_DS1[[j]][[k]]$Symbol %in% "NA")))
    print(table(is.na(data_DS1[[j]][[k]]$Symbol)))
    print("After:")
    data_DS1[[j]][[k]]$Symbol=getSymbols(data_DS1[[j]][[k]]$EnsID, ens)
    print(length(which(data_DS1[[j]][[k]]$Symbol %in% "NA")))
    print(table(is.na(data_DS1[[j]][[k]]$Symbol)))
    print("_____________________________")
 
  }
for (j in c("circData_DS2"  , "sjData_DS2"   ,  "geneData_DS2"))
  for (k in c(1: length(data_DS2[[j]])))
  {  print(paste(j, "____", k))
    print("Before:")
    print(length(which(data_DS2[[j]][[k]]$Symbol %in% "NA")))
    print(table(is.na(data_DS2[[j]][[k]]$Symbol)))
    print("After:")
    data_DS2[[j]][[k]]$Symbol=getSymbols(data_DS2[[j]][[k]]$EnsID, ens)
    print(length(which(data_DS2[[j]][[k]]$Symbol %in% "NA")))
    print(table(is.na(data_DS2[[j]][[k]]$Symbol)))
    print("_____________________________")
   
  }
#################################################################### Save
save(data_DS1, file="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/DATA_TABLES/data_DS1.rda")
save(data_DS2, file="/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/DATA_TABLES/data_DS2.rda")
write.csv(sampleInfo_DS1, "/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/DATA_TABLES/sampleInfo_DS1.csv", row.names=FALSE)
write.csv(sampleInfo_DS2, "/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/DATA_TABLES/sampleInfo_DS2.csv", row.names=FALSE)




