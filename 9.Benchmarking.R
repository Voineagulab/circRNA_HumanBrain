rm(list=ls())
###########
makeSTARsummary=function(path) 
{
  # Makes a STAR summary file based on a directory with standard STAR output
  
  setwd(path)
  files <- dir(pattern = "^.*Log.final.out$" , recursive=TRUE)
  files=files[-grep("STARpass1", files)]
  files=files[-grep("mate", files)]
  logs <- vector(mode = "list")
  
  for (i in 1:length(files))
  {
    logs[[i]] <- read.table(files[i], as.is = TRUE, sep = "\t", fill = TRUE)
    names(logs)[i] <- files[i]
    logs[[i]][, 1] <- trimws(logs[[i]][, 1])
  }
  names(logs) <- sub("Log.final.out$", "", names(logs))
  
  summary=logs[[1]][,1]
  for (j in c(1: length(files)))
  {summary=cbind(summary, logs[[j]][,2])
  colnames(summary)[ncol(summary)]=names(logs)[j]
  }
  return(summary)
}
##Get library size from the STAR output
info=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/NegControl_DATA/SampleInfo.csv")
summaryTable=makeSTARsummary("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/NegControl_DATA/RESULTS/STAR_Output/")
write.table(summaryTable, "/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/NegControl_DATA/RESULTS/STAR_Output/STARsummary.txt", sep="\t")
n=grep("Uniquely mapped reads number |", summaryTable[,1], fixed=TRUE)

##Combine DCC output for all samples
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/NegControl_DATA/RESULTS/DCC_RESULTS/DCC_individual_sample_run/")

info$Label=paste(info$BrainID, info$BrainRegion, info$LibraryPrep, sep="_")
info$nCircDetected=NA
info$nCountEntries=NA
info$totalNreads=NA

##Generate a summary table for the DCC data 
for (j in c(1:nrow(info)))
{
  coord=read.table(paste("CircCoordinates_", info$FileName[j], "_final", sep=""), sep="\t", header=TRUE)
  counts=read.table(paste("CircRNACount_", info$FileName[j], "_final", sep=""), sep="\t", header=TRUE)
  
  info$nCircDetected[j]=nrow(coord)
  info$nCountEntries[j]=nrow(counts)
  
  k=grep(info$FileName[j], colnames(summaryTable))
  info$totalNreads[j]=summaryTable[n,k]
  
  if (j==1) circdata=coord else circdata=rbind(circdata, coord)
}
info$totalNreads=as.numeric(info$totalNreads)
info$circ_per_million=info$nCircDetected/(info$totalNreads/10^6)
write.csv(info, "/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/NegControl_DATA/SampleInfo_annotated.csv")

aggregate(info$nCircDetected ,by=as.data.frame(info$LibraryPrep), FUN=mean)
aggregate(info$totalNreads ,by=as.data.frame(info$LibraryPrep), FUN=mean)
aggregate(info$circ_per_million ,by=as.data.frame(info$LibraryPrep), FUN=mean)
wilcox.test(info$totalNreads ~ info$LibraryPrep)

## Generate a circRNA counts table
#Format consistently with DS1
circdata$Chr=paste0("chr", circdata$Chr)
circdata$Coordinate=paste(circdata$Chr, circdata$Start, circdata$End, circdata$Strand, sep="_")
circdata$Id=paste(circdata$Chr, circdata$Start, circdata$End,  sep="_")

keep=unique(circdata$Coordinate)
circdata=circdata[match(keep, circdata$Coordinate) , ]

for (j in c(1:nrow(info)))
{
  coord=read.table(paste("CircCoordinates_", info$FileName[j], "_final", sep=""), sep="\t", header=TRUE)
  coordinate=paste(paste0("chr", coord$Chr), coord$Start, coord$End, coord$Strand, sep="_")
  counts=read.table(paste("CircRNACount_", info$FileName[j], "_final", sep=""), sep="\t", header=TRUE)
  circdata=cbind(circdata, counts[match(circdata$Coordinate , coordinate), ncol(counts)])
  colnames(circdata)[ncol(circdata)]=info$Label[j]
}  

circdata=as.data.frame(circdata)
infocols=c(1:10)
samplecols=c(11:20)
for (j in samplecols)
  circdata[which(is.na(circdata[,j])==TRUE) , j]=0

# filter for 2 reads in 2 samples
circCounts=circdata[rowSums(circdata[, samplecols] >=2 ) >=2 , ]
circCpm=circCounts

# CalculateCPM
libsize=info$totalNreads[match(colnames(circCounts)[samplecols] , info$Label)]/10^6
for(j in c(1:length(samplecols)))
  circCpm[, samplecols[j]]=circCounts[, samplecols[j]]/libsize[j]

circNegCTL=list(sampleInfo=info, libsize=libsize, circCounts=circCounts, circCpm=circCpm)
save(circNegCTL, file="../../circNegCTL.rda")

# Separate the data for PolyA selected and Ribodepleted
circCpm_PA=circCpm[, c(infocols, grep("_PA", colnames(circCpm)))]
circCpm_RD=circCpm[, c(infocols, grep("_RD", colnames(circCpm)))]
names=gsub("_RD", "", colnames(circCpm_RD)[-infocols])

#Compare with DS1 data
load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/DATA_TABLES/data_DS1.rda")

ds1=data_DS1$circData_DS1$circCounts
length(intersect(circCpm$Coordinate, ds1$Coordinate))/nrow(circCpm)
length(intersect(circCpm$Id, ds1$Id))/nrow(circCpm)

for (s in c(1:length(names)))
{
k=grep(names[s], colnames(ds1)); 
if (s==1) matched_samples=k[1] else matched_samples=c(matched_samples, k[1]) 
}

ds1_matched=ds1[,c(1:15, matched_samples) ]
ds1_matched=ds1_matched[rowSums(ds1_matched[, -c(1:15)] >0 ) >=1 , ]
snames=gsub("_PA", "", colnames(circCpm_PA)[-infocols])
colnames(ds1_matched)[-c(1:15)]=colnames(circCpm_RD)[-infocols]=colnames(circCpm_PA)[-infocols]=snames

results=matrix(0, ncol=5, nrow=11)
colnames(results)=snames
rownames(results)=c("nCirc_DS1", "libsize_DS1", 
                    "nCirc_PA", "libsize_PA", 
                    "nCirc_RD", "libsize_RD",
                    "nCirc_intersectDS1.PA", "nCirc_intersectRD.PA","nCirc_intersectDS1.RD.PA",
                    "fdr_DS1","fdr_RD" )
results=as.data.frame(results)

for (j in c(1:length(snames)))
{
results["nCirc_DS1", j]=length(which(ds1_matched[,match(snames[j], colnames(ds1_matched))] > 0))
results["nCirc_PA", j]=length(which(circCpm_PA[,match(snames[j], colnames(circCpm_PA))] > 0))
results["nCirc_RD", j]=length(which(circCpm_RD[,match(snames[j], colnames(circCpm_RD))] > 0))
results["libsize_DS1", j]=data_DS1$libsizes_DS1[grep(snames[j], names(data_DS1$libsizes_DS1))[1]]
results["libsize_PA", j]=info$totalNreads[grep(paste0(snames[j], "_PA"), info$Label)]/10^6
results["libsize_RD",j]=info$totalNreads[grep(paste0(snames[j],"_RD"), info$Label)]/10^6
DS1.PA=intersect(ds1_matched$Coordinate[which(ds1_matched[,match(snames[j], colnames(ds1_matched))] > 0)],
                 circCpm_PA$Coordinate[which(circCpm_PA[,match(snames[j], colnames(circCpm_PA))] >0 )])
results["nCirc_intersectDS1.PA", j]=length(DS1.PA)
RD.PA=intersect(circCpm_RD$Coordinate[which(circCpm_RD[,match(snames[j], colnames(circCpm_RD))] > 0)],
                circCpm_PA$Coordinate[which(circCpm_PA[,match(snames[j], colnames(circCpm_PA))] >0 )])
results["nCirc_intersectRD.PA", j]=length(RD.PA)
results["nCirc_intersectDS1.RD.PA", j]=length(intersect(DS1.PA, RD.PA))
}

results["fdr_DS1", ]=(results["nCirc_intersectDS1.PA", ]/results["libsize_PA", ])/(results["nCirc_DS1", ]/results["libsize_DS1",])
results["fdr_RD", ]=(results["nCirc_intersectRD.PA", ]/results["libsize_PA", ])/(results["nCirc_RD", ]/results["libsize_RD",])
write.csv(results, "../../nCirc.csv")
# Supplementary Figure 2 B-D is generated by plotting the data from nCirc.csv in Excel

##How many circRNAs have ambiguous strand?
##DS1:
dim(data_DS1$circData_DS1$circCpm_filter)
#[1] 14386   159
length(grep("_.", data_DS1$circData_DS1$circCpm_filter$Coordinate, fixed=TRUE))
#182

load("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/DATA_TABLES/data_DS2.rda")
dim(data_DS2$circData_DS2$circCpm_filter)
#[1] 9440   68
length(grep("_.", data_DS2$circData_DS2$circCpm_filter$Coordinate, fixed=TRUE))
#[1] 54



