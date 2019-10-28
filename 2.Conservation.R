## NOTE: No changes for revision
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
rm(list=ls())

############################################ LOAD DATA ############################################
load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
circCts_ds1=data_DS1$circData_DS1$circCounts_filter; circCpm_ds1=data_DS1$circData_DS1$circCpm_filter
circCts_ds2=data_DS2$circData_DS2$circCounts_filter; circCpm_ds2=data_DS2$circData_DS2$circCpm_filter

bed_ds1=circCts_ds1[, c("Chr", "Start", "End", "Id", "Symbol", "Strand")]
bed_ds2=circCts_ds2[, c("Chr", "Start", "End", "Id", "Symbol", "Strand")]

ds1=paste(circCts_ds1[, 3], ":", circCts_ds1[, 4], "-", circCts_ds1[, 5])
ds2=paste(circCts_ds2[, 3], ":", circCts_ds2[, 4], "-", circCts_ds2[, 5])

write.table(bed_ds1, "Analysis/Mouse_vs_Human/circ_filter_ds1.bed",sep="\t",  row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bed_ds2, "Analysis/Mouse_vs_Human/circ_filter_ds2.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Boxplots of mean expression of conserved vs nonconserved circRNAs
# any novel ones among the conserved?

write.table(bed_ds1[c(1:7000), ], "Analysis/Mouse_vs_Human/circ_filter_ds1_part1.bed",sep="\t",  row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bed_ds1[c(7001: 14386), ], "Analysis/Mouse_vs_Human/circ_filter_ds1_part2.bed",sep="\t",  row.names=FALSE, col.names=FALSE, quote=FALSE)
############################################ Liftover in UCSC genome browser

############################################
ds1_1=read.table("Analysis/Mouse_vs_Human/circ_filter_ds1_part1_liftover_mm9.bed", sep="\t")
ds1_2=read.table("Analysis/Mouse_vs_Human/circ_filter_ds1_part2_liftover_mm9.bed", sep="\t")
ds2=read.table("Analysis/Mouse_vs_Human/circ_filter_ds2_liftover_mm9.bed", sep="\t")
mm9_circ=read.csv("Analysis/Mouse_vs_Human/Rybak-Wolf_suppTables_mouseCircRNA.csv")

ds1=rbind(ds1_1, ds1_2)

colnames(ds1)=c("mm9_Chr", "mm9_Start", "mm9_End", "hg19_Id", "Symbol", "Strand")
colnames(ds2)=c("mm9_Chr", "mm9_Start", "mm9_End", "hg19_Id", "Symbol", "Strand")

ds1$mm9_Id=paste(ds1$mm9_Chr, ds1$mm9_Start, ds1$mm9_End, sep="_")
ds1$mm9_Coordinate=paste(ds1$mm9_Chr, ds1$mm9_Start, ds1$mm9_End, ds1$Strand, sep="_")

ds2$mm9_Id=paste(ds2$mm9_Chr, ds2$mm9_Start, ds2$mm9_End, sep="_")
ds2$mm9_Coordinate=paste(ds2$mm9_Chr, ds2$mm9_Start, ds2$mm9_End, ds2$Strand, sep="_")

mm9_circ=read.csv("Analysis/Mouse_vs_Human/Rybak-Wolf_suppTables_mouseCircRNA.csv")
mm9_circ$start=mm9_circ$start+1 ; 
mm9_circ$Id=paste(mm9_circ[,1],mm9_circ[,2],mm9_circ[,3], sep="_")
mm9_circ$Coordinate=paste(mm9_circ[,1],mm9_circ[,2],mm9_circ[,3], mm9_circ[,6] ,sep="_")

length(intersect(ds1$mm9_Id, mm9_circ$Id))
#[1] 1195
length(intersect(ds1$mm9_Coordinate, mm9_circ$Coordinate))
#[1] 1190
length(intersect(ds2$mm9_Coordinate, mm9_circ$Coordinate))
#[1] 880
length(intersect(ds2$mm9_Id, mm9_circ$Id))
#[1] 883
length(intersect(intersect(ds1$mm9_Coordinate, mm9_circ$Coordinate), intersect(ds2$mm9_Coordinate, mm9_circ$Coordinate)))
#[1] 834
#Rybak:Mouse:15849; Human:65731; Conserved:4523 --> ~ 6% of the human ones are conserved
#Our data: 8.6% in  DS11; 9.6%in DS2

length(intersect(intersect(ds1$mm9_Coordinate, mm9_circ$Coordinate), intersect(ds2$mm9_Coordinate, mm9_circ$Coordinate)))
# 94% overlap between the conserved ones between ds1 and ds2

######### Vs You et al. NN 2015
mm9_circ=read.csv("Analysis/Mouse_vs_Human/You_NN2015_mouse_circRNA_TPM.csv")
mm9_circ$Gname=paste0("chr", mm9_circ$Gname)
temp=transpose(strsplit(mm9_circ$Gname, "_", fixed=TRUE))
## So frustrating!!! they have the end coord before the start coord! wasted 1h.
mm9_circ$Id=paste(temp[[1]], as.numeric(temp[[3]]), as.numeric(temp[[2]]), sep="_")
length(intersect(mm9_circ$Id, ds1$mm9_Id))
mm9_circ2=read.csv("Analysis/Mouse_vs_Human/Rybak-Wolf_suppTables_mouseCircRNA.csv")
mm9_circ2$start=mm9_circ2$start+1 ; 
mm9_circ2$Id=paste(mm9_circ2[,1],mm9_circ2[,2],mm9_circ2[,3], sep="_")
length(intersect( mm9_circ2$Id, ds1$mm9_Id))
       
length(intersect(ds1$mm9_Id, mm9_circ$Id))
#[1] 1866 
length(intersect(ds2$mm9_Id, mm9_circ$Id))
#[1] 1304
dim(mm9_circ)
#[1] 28482    31
#> 1866/nrow(ds1)
#[1] 0.1349631
#> 1304/nrow(ds2)
#[1] 0.1423736
       