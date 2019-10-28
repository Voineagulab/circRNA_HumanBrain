rm(list=ls())
library(WGCNA)
library(multtest)
library(gplots)
library(goseq)
options(stringsAsFactors = FALSE)
source("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/SCRIPTS/Final/Supplementary_Code/Functions.R")
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
############# LOAD DATA 
load("DATA_TABLES/data_DS1.rda")
load("DATA_TABLES/ratios_ds1.rda")

samples1=data_DS1$sampleInfo_DS1
libsizes1=data_DS1$libsizes_DS1
infocolsC=c(1:15)

##########Select CTX samples
CTX=as.character(samples1$Sample[which(samples1$BroadRegionID=="CTX")])
samples1=samples1[CTX,]; libsizes1=libsizes1[CTX]

##########Select circRNAs expressed in at least half of the samples
circ=ratios_ds1$ci[, CTX]
keep=which(rowSums(circ[,] > 0) >= ncol(circ)/2)
circ=circ[ keep, ]

##########Calculate residuals on log2-transformed data to remove the effect of covariates
datCTX=residuals(lm(t(log2(circ+0.5)) ~ Neurons + Sex +  SeqBatch + BrainBank + Age + Median.3prime.Bias.picard, samples1 ))

##########Select soft-thresholding beta power
powers=c(c(1:10), seq(from = 12, to = 20, by = 2))  #This is the list of powers we will be testing
sft = pickSoftThreshold(datCTX, powerVector = powers, verbose = 5, blockSize = 10000, networkType = "signed")
pdf(file="Analysis/WGCNA/WGCNA_CTX_SoftThreshold.pdf")
 cex1 = 0.9
 par(mfrow = c(1,2))
  #Fit to scale-free topology
 plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",
      ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"));
 text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red");
 abline(h=c(0.9 , 0.8, 0.5), col="red")
  #Mean connectivity
 plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold (power)",
      ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity"))
 text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col="red")
dev.off()

##########Construct the network using beta-power = 10
net=blockwiseModules(datCTX, power=10, numericLabels=TRUE, networkType = "signed", corFnc="bicor", 
                     minModuleSize=10, mergeCutHeight=0.2, saveTOMs=FALSE, verbose=6, 
                     nThreads=24, maxBlockSize=35000, checkMissingData=FALSE)
adjacency = adjacency(datCTX, corFnc = "bicor", type = "signed", power = 10)
rownames(adjacency)=colnames(adjacency)=colnames(datCTX)

modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N") 
modules$Label=paste("M", modules$Label, sep=""); 
modules$Color=c("grey",labels2colors(modules$Label[-1])) 
moduleLabel=paste("M",net$colors, sep=""); 
moduleColor=modules$Color[match(moduleLabel, modules$Label)] 

########## Calculate kMEs
KMEs<-signedKME(datCTX, net$MEs,outputColumnName = "M") 
rownames(KMEs)=colnames(datCTX)
kme=data.frame("NA", moduleColor,moduleLabel, KMEs)
colnames(kme)[1]="Symbol"
c=grep("chr", rownames(kme)); kme$Symbol[c]=as.character(data_DS1$circData_DS1$circCpm$Symbol[match(rownames(kme)[c], data_DS1$circData_DS1$circCpm$Coordinate)])

########## Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.1, will be assigned to the grey/M0 module
kmeInfoCols=c(1:3)
kmedata=kme[,-kmeInfoCols];
pvalBH=kmedata; pvalBH[,]=NA
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedata[,j], nSamples=ncol(datCTX)), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}
kme$newModule="NA"
for (j in c(1:nrow(kmedata)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m=which(kmedata[j,]==max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedata[j,m]>0.1)) kme$newModule[j]=as.character(colnames(kmedata)[m])
}
#Assign genes not associated to any module to M0
kme$newModule[which(kme$newModule%in%"NA")]="M0"
#Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)]
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor
kme=kme[,-grep("newModule", colnames(kme))];kme=kme[,-grep("newColor", colnames(kme))]

##########Save kMEs 
mod=modules$Label[-1] 
kmeTable=kme[,kmeInfoCols]
for(j in c(1:length(mod)))
{
  kmeTable=cbind(kmeTable, kmedata[,match(mod[j],colnames(kmedata))]);colnames(kmeTable)[ncol(kmeTable)]=paste("kME", mod[j], sep="_")                                                                             
  kmeTable=cbind(kmeTable, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTable)[ncol(kmeTable)]=paste("pvalBH", mod[j], sep="_")
}
write.csv(kmeTable, "Analysis/WGCNA/kmE_CTX.csv")

##########Save Module Eigengenes 
me<-data.frame(rownames(datCTX), net$MEs) # Sample names bound to module eigengenes
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me)[1]="Sample"
write.csv(me, "Analysis/WGCNA/ME_CTX.csv", row.names=FALSE)

########## Plot Module Eigengene Values
age=samples1$Age[match(me$Sample, samples1$Sample)]
pheno=samples1$ASD.CTL[match(me$Sample, samples1$Sample)]
np=samples1$Neurons[match(me$Sample, samples1$Sample)]
seizures=samples1$Seizures[match(me$Sample, samples1$Sample)]
colors=labels2colors(pheno, colorSeq = c("salmon", "lightblue"))

pdf("Analysis/WGCNA/moduleBarplots_CTX.pdf", height=5, width=15)
par(mfrow=c(1,2))
mod=paste("M", c(0:(ncol(me)-2)), sep="")
genes=grep("ENS", rownames(kme))
circ=grep("chr", rownames(kme))
########## Statistical tests
stats=as.data.frame(matrix(NA, nrow=length(mod), ncol=6))
rownames(stats)=mod
colnames(stats)=c("corAge", "corAge.p", "corNP", "corNP.p", "wilcox.test.ASD", "wilcox.test.Seizures")
for(m in c(1:length(mod))) 
{ 
  j=match(mod[m], colnames(me))
  k=grep(mod[m], kme$moduleLabel)
  stats$corAge[m]=cor(me[,j], age, method="s")
  stats$corAge.p[m]=corPvalueStudent(stats$corAge[m], nrow(datCTX)) * (nrow(stats) - 1)
  stats$corNP[m]=cor(me[,j], np, method="s")
  stats$corNP.p[m]=corPvalueStudent(stats$corNP[m], nrow(datCTX))* (nrow(stats) - 1)
  stats$wilcox.test.ASD[m]=wilcox.test(me[,j]~pheno)$p.value * (nrow(stats) - 1)
  stats$wilcox.test.Seizures[m]=wilcox.test(me[,j]~seizures)$p.value * (nrow(stats) - 1)
}
write.csv(stats, "Analysis/WGCNA/stats.csv")
for(m in mod) 
{ 
  j=match(m, colnames(me)); k=grep(m, kme$moduleLabel); s=grep(m, rownames(stats))
  title=paste(m, ", Total=", length(k), ", ncirc=", length(intersect(circ,k)))
  text=paste("corAge=", stats$corAge[s], ",corAge.p=", stats$corAge.p[s], "corNP=", stats$corNP[s], ",corNP.p=", stats$corNP.p[s],", wilcoxASD.p=", stats$w[s])
  order=c(which(pheno=="ASD"), which(pheno=="CTL"))
  barplot(me[order,j],  xlab="Samples", ylab="ME",col=colors[order], main=title, names="", axisnames = FALSE, sub=text)
  boxplot(me[,j]~ pheno,   ylab="ME")
}
dev.off()

########## Export data to Cytoscape for M4: top 20 nodes;top 50 edges.
m4=rownames(kme)[which(kme$moduleLabel=="M4")]
kme_m4=kme[m4,]
kme_m4_sorted=kme_m4[order(kme_m4$M4, decreasing=TRUE) ,]

adj=adjacency[rownames(kme_m4_sorted)[1:20],rownames(kme_m4_sorted)[1:20]]
rownames(adj)=colnames(adj)=paste(rownames(kme_m4_sorted)[1:20], kme_m4_sorted$Symbol[1:20], sep="_")
cytoscape=exportNetworkToCytoscape(adj, weighted = TRUE ,  threshold = 0) 

edges=cytoscape$edgeData  
edgesSorted=edges[order(edges$weight, decreasing=TRUE) ,]
write.table(edgesSorted[1:50, ], "Analysis/WGCNA/m4_topEdges.txt", sep="\t", row.names=FALSE, quote=FALSE)

########## Module Eigengene Boxplot for M4
pdf("Analysis/WGCNA/Fig7.m4.pdf", height=4, width=8)
par(mfrow=c(1,3))
boxplot(me$M4 ~ pheno, col=c("dodgerblue2", "gold2"), notch=TRUE, ylab="Module Eigengene Value")
dev.off()

############

# load("DATA_TABLES/data_DS1.rda")
# genes=data_DS1$geneData_DS1$rpkm_filter
# me=read.csv("Analysis/WGCNA/ME_CTX.csv")
# 
# rownames(me)=me[,1]; me=me[,-1]
# genes_CTX=t(log2(genes[, match(rownames(me), colnames(genes))])+0.5)
# 
# cormx=t(cor(me, genes_CTX, method="s"))
# cormx=data.frame(as.character(genes$Symbol[match(rownames(cormx), rownames(genes))]), cormx)
# colnames(cormx)[1]="Symbol"


#####Module Preservation

load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/ratios_ds2.rda")
samples2=data_DS2$sampleInfo_DS2
libsizes2=data_DS2$libsizes_DS2

CTX2=as.character(samples2$Sample[which(samples2$BroadRegionID=="CTX")])
samples2=samples2[CTX2,]; libsizes2=libsizes2[CTX2]
circ2=ratios_ds2$ci[, CTX2]
keep2=which(rowSums(circ2[,] > 0) >= ncol(circ2)/2)
circ2=circ2[ keep2, ]
datCTX2=residuals(lm(t(log2(circ2+0.5)) ~ Neurons + Sex + RIN  + Age + Median.3prime.Bias.picard, samples2 ))

multiExpr=list(DS1=list(data=datCTX), DS2=list(data=datCTX2))
labels=kmeTable$moduleColor; names(labels)=rownames(kmeTable)
multiColor=list(DS1=t(labels))
mp = modulePreservation(multiExpr, multiColor, dataIsExpr = TRUE, networkType = "signed",checkData = FALSE,
                        referenceNetworks = 1,
                        nPermutations = 200,
                        randomSeed = 1,
                        verbose = 3)
save(mp, file = "Analysis/WGCNA/modulePreservation.RData");
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot

pdf("Analysis/WGCNA/ModulePreservation.pdf", width=10, height=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()
###Enrichment analysis
m4Genes=as.character(kmeTable$Symbol[which(kmeTable$moduleLabel%in%"M4")])
expGenes=as.character(data_DS1$geneData_DS1$rpkm_filter$Symbol)
###Gprofiler

library(gprofiler2)
go_path=gost(m4Genes, organism = "hsapiens", ordered_query = FALSE,
             multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
             measure_underrepresentation = FALSE, evcodes = FALSE,
             user_threshold = 0.05, correction_method ="fdr", custom_bg = expGenes,
             numeric_ns = "")
m4_go=go_path$result
write.csv(m4_go, "Analysis/WGCNA/m4_go.csv")
# Cell type Markers
mk=read.csv("Annotation_GeneLists/PEC_DER-19_Single_cell_markergenes_TPM.csv")
mkMod=unique(mk$CellType)
m4Mod=kmeTable$Symbol[which(kmeTable$moduleLabel%in%"M4")]
p.hyper=rep(0, length(mkMod)); 
names(p.hyper)=mkMod
for (j in c(1:length(mkMod)))
{
  mkM=mk$GeneName[which(mk$CellType%in%mkMod[j])]
  p.hyper[j]=phyper(length(intersect(mkM, m4Mod)),
                    length(m4Mod),
                    (length(expGenes)-length(m4Mod)),
                    length(mkM),
                    lower.tail = FALSE)
}
write.csv(p.hyper, "Analysis/WGCNA/m4_cellTypeMarkers.csv")
#Parikshak modules
pk=read.csv("Annotation_GeneLists/Parikshak_CTX_SuppTable.csv")
pkMod=unique(pk$WGCNA.Module.Label)
m4Mod=intersect(pk$HGNC.Symbol, kmeTable$Symbol[which(kmeTable$moduleLabel%in%"M4")])
common=intersect(pk$HGNC.Symbol, expGenes)
  p.hyper=rep(0, length(pkMod)); 
  names(p.hyper)=pkMod
  for (j in c(1:length(pkMod)))
  {
    pkM=intersect(pk$HGNC.Symbol[which(pk$WGCNA.Module.Label%in%pkMod[j])],common)
    
    p.hyper[j]=phyper(length(intersect(pkM, m4Mod)),
                      length(m4Mod),
                      (length(common)-length(m4Mod)),
                      length(pkM),
                      lower.tail = FALSE)
  }
write.csv(p.hyper, "Analysis/WGCNA/m4_ParikshakModuleOverlap.csv")
#No overlap with Parikshak et al modules that were differentially expressed in ASD
rm(list=ls())
kme=read.csv("Analysis/WGCNA/kmE_CTX.csv")
load("DATA_TABLES/data_DS1.rda")

m4=which(kme$moduleLabel%in%"M4")
circGenes_Symbol=unique(kme$Symbol[m4])
geneR=unique(data_DS1$geneData_DS1$rpkm_filter$Symbol)

calculate.overlaps=function(files)
{
  #each file in the list should contain a gene category. The function calculates the significance of the overlap with circRNA expressing genes.
  #hyper-geometric test p-values are Bonferroni corrected.
  p.hyper=rep(0, length(files)); prop=rep(0, length(files)); prop.exp=rep(0, length(files)); mg=rep(" ", length(files))
  names(p.hyper)=names(prop)=names(prop.exp)=gsub("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/", "", files, fixed=TRUE)
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
ndd_files=paste0("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/Annotation_GeneLists/", c("ID_Vissers_et_al. 2016_supp_table.txt", "SFARI-15-06-2018_Syndromic_or_Score_1_4.txt", "SZDB_minScore2.txt"))
ndd_overlaps=calculate.overlaps(ndd_files)
write.csv(ndd_overlaps, "Analysis/WGCNA/M4_NDD_overlaps.csv")
