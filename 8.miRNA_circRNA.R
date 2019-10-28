library(data.table)
library(ggplot2)
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
load("../DATA_TABLES/data_DS1.rda")
annot=read.csv("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Analysis/DatasetCharacterisationPlots/SuppTable_DS1_annot.csv")
mirna=read.csv("../../PublishedDatasets/Wu_nn.4373-S4_miRNARawCounts.csv")
rownames(mirna)=mirna$miR.ID
mirna=mirna[,-1]
ntotal=apply(mirna,2,sum)
mirnaN=mirna; mirnaN[,]=NA
for(j in c(1:ncol(mirnaN))) mirnaN[,j]=mirna[,j]/(ntotal[j]/10^6)
rm(mirna)

circ=data_DS1$circData_DS1$circCpm_filter
infocolsC=c(1:15)
samples1=data_DS1$sampleInfo_DS1
# Select common samples (n=95)
m=match(colnames(circ), colnames(mirnaN))
keep=which(is.na(m)==FALSE)
circUse=circ[,keep]
circUse=circUse[rowSums(circUse > 0 ) > ncol(circUse)/2,]
mirnaUse=mirnaN[,m[keep]]
genesUse=data_DS1$geneData_DS1$rpkm_filter[, colnames(circUse)]

circUse=residuals(lm(t(log2(circUse +0.05)) ~  Sex + Age + Median.3prime.Bias.picard + SeqBatch + BrainBank, samples1[colnames(circUse), ]   ))
mirnaUse=residuals(lm(t(log2(mirnaUse +0.05)) ~   Sex + Age + Median.3prime.Bias.picard + SeqBatch + BrainBank, samples1[colnames(mirnaUse), ]  ))
genesUse=residuals(lm(t(log2(genesUse +0.5)) ~   Sex + Age + Median.3prime.Bias.picard + SeqBatch + BrainBank, samples1[colnames(genesUse), ]  ))
cyrano=grep("ENSG00000247556", colnames(genesUse))
cdr=grep("chrX_139865340_139866824_.", colnames(circUse))
mir7=grep("hsa-miR-7-5p", colnames(mirnaUse))
mir671=grep("hsa-miR-671", colnames(mirnaUse))
plotdata=cbind(genesUse[,cyrano], circUse[,cdr] , mirnaUse[,mir7], mirnaUse[,mir671[1]] ,mirnaUse[,mir671[2]], samples1[rownames(circUse),])

# cyrano=grep("ENSG00000247556", rownames(genesUse))
# cdr=grep("chrX_139865340_139866824_.", rownames(circUse))
# mir7=grep("hsa-miR-7-5p", rownames(mirnaUse))
# mir671=grep("hsa-miR-671", rownames(mirnaUse))
# plotdata=cbind(t(log10(genesUse[cyrano, ])), t(log10(circUse[cdr, ])) , t(log10(mirnaUse[mir7,])), t(log10(mirnaUse[mir671[1], ])) ,t(log10(mirnaUse[mir671[2], ])) , samples1[colnames(circUse),])

plotdata=as.data.frame(plotdata)
colnames(plotdata)[1:5]=c("cyrano", "CDR1AS", "mir7_5p", "mir671_3p", "mir671_5p")

pdf("Analysis/miRNA_circRNA.pdf", height=4,width=4)
ggplot(plotdata, aes(x=CDR1AS, y=cyrano, colour=mir7_5p)) +geom_point()  + geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  scale_colour_gradient2(high = "red", mid = "white",
                         low  = "blue", midpoint = mean(plotdata$mir7_5p), space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")
ggplot(plotdata, aes(x=CDR1AS, y=mir671_5p, colour=mir7_5p)) +geom_point() + geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  scale_colour_gradient2(high = "red", mid = "white",
                         low  = "blue", midpoint = mean(plotdata$mir7_5p), space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour")
dev.off()
cor(plotdata$cyrano,plotdata$CDR1AS, method="s" )
#[1,]0.6368281
cor(plotdata$mir671_5p,plotdata$CDR1AS, method="s" )
#[1] -0.5765258
cor(plotdata$mir7_5p,plotdata$CDR1AS, method="s" )
#[1] -0.6368421
cor(plotdata$cyrano,plotdata$mir7_5p, method="s" )
#[1] -0.6677352

circ_mirna=cor(circUse, mirnaUse, method="s", use="pair")
circ_mirna_id=annot$circBase_circRNA.ID[match(rownames(circ_mirna), annot$Coordinate)]

annot=read.csv("Analysis/DatasetCharacterisationPlots/SuppTable_DS1_annot.csv")
starbase=read.table("Annotation_GeneLists/starBaseV3_hg19_CLIP-seq_all_circ.txt", sep="\t", header=TRUE)
starbase=starbase[,c(1:5,10)]
circ_mirna_sb=circ_mirna[which(circ_mirna_id%in%starbase$geneID), 
                         which(colnames(circ_mirna)%in%starbase$miRNAname)]

rownames(circ_mirna_sb)=annot$circBase_circRNA.ID[match(rownames(circ_mirna_sb), annot$Coordinate)]

interactions=as.data.frame(melt(circ_mirna_sb))
colnames(interactions)=c("circRNA", "miRNA", "SpearmanCor")
interactions$intID=paste(interactions$circRNA, interactions$miRNA, sep="_")

starbase=starbase[which((starbase$geneID%in%circ_mirna_id)&(starbase$miRNAname%in%colnames(circ_mirna))), ]
starbase$intID=paste(starbase$geneID, starbase$miRNAname, sep="_")
length(unique(starbase$miRNAname))
#[1] 282
length(unique(starbase$geneID))
#[1] 147
m=match(starbase$intID, interactions$intID)
starbase$SpearmanCor=interactions$SpearmanCor[m]

keep=which(abs(starbase$SpearmanCor)> 0.5)
starbase_th=starbase[keep,]
length(unique(starbase_th$miRNAname))
#[1] 63
length(unique(starbase_th$geneID))
#[1] 31
starbase_th$geneSymbol=annot$Symbol[match(starbase_th$geneID, annot$circBase_circRNA.ID)]
starbase_th$Coord=annot$Coordinate[match(starbase_th$geneID, annot$circBase_circRNA.ID)]
write.csv(starbase_th, "Analysis/mirna_circrna_interactions.csv")
length(unique(starbase_th$intID))
#101
m=match(unique(starbase_th$intID),starbase_th$intID)
table(sign(starbase_th$SpearmanCor[m]))
#-1  1 
#64 37
table=starbase_th[m,]
write.csv(table, "Analysis/mirna_circrna_interactions_SupTable.csv")
