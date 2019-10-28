
############################################  PLOTTING FUNCTIONS
######################################################################################################################################################
plotVennCirc=function (DS1, DS2, circBase, title)
{library(VennDiagram)
  area1=nrow(DS1); area2=nrow(DS2); area3=nrow(circBase)
  n12=length(intersect(DS1$Coordinate, DS2$Coordinate))
  n23=length(intersect(DS2$Coordinate, circBase$Coordinate))
  n13=length(intersect(DS1$Coordinate, circBase$Coordinate))
  n123=length(intersect(intersect(DS1$Coordinate, DS2$Coordinate), circBase$Coordinate))
  venn.plot <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, alpha=c(0.5, 0.5, 0.5), fill=c("mediumorchid", "seagreen3", "tan1"), category=c("DS1", "DS2", "circBase"), main=title)
  grid.draw(venn.plot);grid.newpage()
  venn.plot <- draw.pairwise.venn(area1, area2,n12, euler.d=TRUE, scaled = TRUE,  alpha=c(0.5, 0.5),fill=c("mediumorchid", "seagreen3"), category=c("DS1", "DS2"), main=title)
  grid.draw(venn.plot);grid.newpage()
}
######################################################################################################################################################
plot4wayVennCirc=function (DS1, DS2, title)
{library(VennDiagram)
  circBase_brain=circBase[grep("Rybak", circBase$circRNA.study),]
  area1=nrow(DS1); area2=nrow(DS2); area3=nrow(circBase); area4=nrow(circBase_brain)
  n12=length(intersect(DS1$Coordinate, DS2$Coordinate))
  n23=length(intersect(DS2$Coordinate, circBase$Coordinate))
  n13=length(intersect(DS1$Coordinate, circBase$Coordinate))
  
  n14=length(intersect(DS1$Coordinate, circBase_brain$Coordinate))
  n24=length(intersect(DS2$Coordinate, circBase_brain$Coordinate))
  n34=length(intersect(circBase$Coordinate, circBase_brain$Coordinate))
  
  n123=length(intersect(intersect(DS1$Coordinate, DS2$Coordinate), circBase$Coordinate))
  n124=length(intersect(intersect(DS1$Coordinate, DS2$Coordinate), circBase_brain$Coordinate))
  n134=length(intersect(intersect(DS1$Coordinate, circBase$Coordinate), circBase_brain$Coordinate))
  n234=length(intersect(intersect(DS2$Coordinate, circBase$Coordinate), circBase_brain$Coordinate))
  
  n1234=length(intersect(intersect(intersect(DS1$Coordinate, DS2$Coordinate), circBase$Coordinate), circBase_brain$Coordinate))
  venn.plot <- draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,n34, n123, n124, n134, n234, n1234,
                              euler.d=TRUE, scaled = TRUE,  alpha=c(0.5, 0.5, 0.5,0.5),
                              fill=c("mediumorchid", "seagreen3", "skyblue", "peachpuff"), 
                              category=c("DS1", "DS2", "circBase_all", "circBase_brain"), 
                              main=title)
  grid.draw(venn.plot);grid.newpage()
}
######################################################################################################################################################
plotPCA=function(data, infocols, samples,title, obs)
{
  samples=samples[match(colnames(data)[-infocols], samples$Sample) , ]
  colors <- factor(samples$RegionID, levels=(c("ba9", "ba41-42-22", "vermis")) ,labels=c("#66c2a5",  "#8da0cb" ,"#fc8d62"))
  p=rep(17, nrow(samples))
  p[which(samples$ASD.CTL %in%"ASD")]=16
  cormx=cor(data[,-infocols], method="s", use="complete")
  pca=princomp(cormx)
  var=round(pca$sdev^2,4)
  plot(pca$loadings[,1], pca$loadings[,2], pch=p, col=as.character(colors), 
       xlab=paste("PC1:", var[1]*100, "%"),
       ylab=paste("PC2:", var[2]*100, "% "), main=title)
  #text(pca[,1], pca[,2], samples$RegionID)
  #legend("topleft", pch=16, col=rainbow_hcl(3), legend=unique(samples$RegionID))
}
######################################################################################################################################################
scatterplot.mean.f=function(data1, data2, xlab, ylab, title, matchcol)
{
  data1=data1[which(data1[, matchcol]%in%data2[,matchcol]) ,]
  m=match(data1[, matchcol], data2[,matchcol])
  print(length(which(is.na(m)==TRUE)))
  corS=cor(data1$Mean, data2$Mean[m],  use="pair", method="s")
  plot(log2(data1$Mean) , log2(data2$Mean[m]),pch=20, col="grey", sub=paste0("Spearman Cor= ",round(corS,2)), xlab=xlab, ylab=ylab, main=title)
  abline(lm(log2(data2$Mean[m]) ~ log2(data1$Mean)), col="red", untf=F)
}
######################################################################################################################################################
heatmap.2.f=function(data,main, symm)
{
  library(gplots)
  heatmap.2(data, main=main,symm=symm,key=TRUE,  hline=0, vline=0, trace="none",tracecol=0, scale="none", margins=c(10,10))
}

############################################  FUNCTIONS FOR LM
######################################################################################################################################################
format_for_lm=function(data, libs, samples, voom)
{ library(limma)
  expdata=data[, match(rownames(samples), colnames(data))]
  #featuredata=data[, -match(rownames(samples), colnames(data))]
  if (voom==TRUE) data_T=voom(expdata, lib.size=libs)$E else data_T=log2(expdata+0.5) 
  data_T=t(data_T)
  return(data_T)
}
######################################################################################################################################################
adjPvals.from.lm=function (coef.dat)
{
  library(data.table)
  coef.dat.t=as.data.frame(transpose(coef.dat))
  nr=nrow(coef.dat[[1]]);  nc=ncol(coef.dat[[1]])
  for(j in c(1:nc)) for (k in c(1:nr)) 
    colnames(coef.dat.t)[(j-1)*nr+k]=paste(colnames(coef.dat[[1]])[j], rownames(coef.dat[[1]])[k], sep="_")
  rownames(coef.dat.t)=names(coef.dat)
  pvals=coef.dat.t[, grep("Pr(>|t|)", colnames(coef.dat.t))]
  # get BH corrected pvals
  pvals.adj=pvals; pvals.adj[,]=NA
  for (j in c(1:ncol(pvals.adj))) pvals.adj[,j]=p.adjust(pvals[,j], method="BH")
  rownames(pvals.adj)=gsub("Response", "", rownames(pvals.adj)); rownames(pvals.adj)=gsub(" ", "", rownames(pvals.adj))
  return(pvals.adj)
}
######################################################################################################################################################
run_lm=function(data,samples,design)
{ 
  model=lm(paste("data ~", design, sep=""),  samples)
  residuals=t(residuals(model))
  adj.pvals=adjPvals.from.lm(coefficients(summary(model)))
  return(list(residuals=residuals, adj.pvals=adj.pvals))
}

############################################  OTHER FUNCTIONS
summaryToLibsizes=function(summary)
{rownames(summary)=summary[,1]
summary=summary[,-1]
colnames(summary)=formatLabels(colnames(summary))
libsizes=t(summary[which(rownames(summary)%in%"Uniquely mapped reads number |"),])
libsizes=as.numeric(libsizes)/10^6
names(libsizes)=colnames(summary)
return(libsizes)
}
######################################################################################################################################################
formatLabels=function(x)
{
  x=gsub("Sample_", "", x)
  x=gsub(".", "_", x, fixed=TRUE)
  x=gsub("-", "_", x, fixed=TRUE)
}
######################################################################################################################################################
reorder.f=function(data.list, sampleinfo, infocols)
{
  for (j in c(1:length(data.list))) 
  { index=match(rownames(sampleinfo), colnames(data.list[[j]]))
  data.list[[j]]=data.list[[j]][,c(infocols,index)]}
  return(data.list)
}
######################################################################################################################################################
calcSummaryData=function(data, infocols, samples, th)
{
  summary=data[, infocols] 
  s=which(colnames(data)%in%samples$Sample)
  data=data[,s]
  summary$Mean=apply(data,1,mean, na.rm=TRUE) 
  summary$Sd=apply(data,1,sd, na.rm=TRUE) 
  summary$Var=apply(data,1,var, na.rm=TRUE) 
  summary$CV=summary$Sd/summary$Mean 
  summary$nexpPerFeature=rowSums(data[,] >= th, na.rm=TRUE)
  summary$nexpPerFeatureCTX=rowSums(data[, grep("ba", colnames(data))] >= th, na.rm=TRUE)
  summary$nexpPerFeatureCB=rowSums(data[,grep("vermis", colnames(data))] >= th, na.rm=TRUE)
  summary$MeanCTX=apply(data[, grep("ba", colnames(data))],1,mean, na.rm=TRUE) 
  summary$MeanVermis=apply(data[, grep("vermis", colnames(data))],1,mean , na.rm=TRUE) 
  
  summary=summariseAnnot(summary)
  return(summary)
}
######################################################################################################################################################
summariseAnnot=function(summary)
{
  summary$summaryAnnot="NA"
  summary$summaryAnnot[grep("intergenic", summary$Annot)]="intergenic"
  summary$summaryAnnot[grep("intronic", summary$Annot)]="intronic"
  summary$summaryAnnot[grep("exonic", summary$Annot)]="exonic"
  summary$summaryAnnot[grep("exonJ", summary$Annot)]="exonJ"
  summary$summaryAnnot[grep("exonJ_exonJ", summary$Annot)]="exonJ-exonJ"
  return(summary)
}

######################################################################################################################################################
runGoSeq=function(geneset, genome="hg19", id="ensGene", lengthbias=1)
{# geneset is a named vector including all genes in the dataset. names: gene ids corresponding to "id"
  # values: 1- DE genes; 0-not DE genes.
  library(goseq)
  pwf <- nullp(geneset, genome, id)
  if (lengthbias){
    GO.wall=goseq(pwf,genome, id,method="Wallenius")
  } else {
    GO.wall=goseq(pwf,genome, id,method="Hypergeometric")
  }
  go <- cbind(GO.wall, p.adjust(GO.wall$over_represented_pvalue, method="BH"))
  colnames(go)[ncol(go)] <- "FDR"
  return(go)
}

######################################################################################################################################################
generate.fasta_for_MEMEChIP=function(input_data, w, out_bed_file, out_fasta_file)
{
  input_data_startBed=as.data.frame(transpose(strsplit(input_data, "_")))
  colnames(input_data_startBed)=c("Chr", "Start", "End", "Strand")
  input_data_startBed$Start=as.numeric(as.character(input_data_startBed$Start)) - w
  input_data_startBed$End=as.numeric(as.character(input_data_startBed$Start))+ 2*w
  input_data_startBed=cbind(input_data_startBed[, c(1:3)],paste0("Start_",input_data), 0, input_data_startBed$Strand )
  
  input_data_endBed=as.data.frame(transpose(strsplit(input_data, "_")))
  colnames(input_data_endBed)=c("Chr", "Start", "End", "Strand")
  input_data_endBed$End=as.numeric(as.character(input_data_endBed$End))+ w
  input_data_endBed$Start=as.numeric(as.character(input_data_endBed$End))-2*w
  input_data_endBed=cbind(input_data_endBed[, c(1:3)],paste0("End_",input_data), 0, input_data_endBed$Strand )
  
  colnames(input_data_startBed)[4:6]=colnames(input_data_endBed)[4:6]=c("Id", "Score", "Strand")
  input_data_bed=rbind(input_data_startBed, input_data_endBed)
  input_data_bed$Chr=gsub("chr", "", input_data_bed$Chr)
  
  input_data_bed$Start=as.integer(input_data_bed$Start)
  input_data_bed$End=as.integer(input_data_bed$End)
  
  write.table(input_data_bed ,out_bed_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  genome=" /Volumes/MacintoshHD_RNA/Users/rna/REFERENCE/HUMAN/Ensembl_GRCh37_hg19/genome/genome.fa"
  runFunction="/Volumes/MacintoshHD_RNA/Users/rna/PROGRAMS/bedtools2/bin/fastaFromBed "
  options="-s -name"
  command=paste0(runFunction, options, " -fi ", genome, " -bed ", out_bed_file, " -fo ", out_fasta_file )
  print(command)
  system(command)
}

# ######################################################################################################################################################
# runDEseq2=function(data, infocols, sample_info,designForm, useSamples, useVar)
# {
#   library(DESeq2)
#   # format sample info matrix
#   sample_info=sample_info[which(rownames(sample_info)%in%useSamples), which(colnames(sample_info)%in%useVar)]
#   #format exp data
#   featuredata=data[, infocols]
#   data=data[, match(rownames(sample_info), colnames(data))]
#   ###Deseq
#   dds <- DESeqDataSetFromMatrix(countData = data,
#                                 colData = sample_info,
#                                 design = designForm)
#   dds <- DESeq(dds)
#   return(dds)
# }
# ######################################################################################################################################################
# getSig=function(res, p, fc)
# {
#   sig=which((res$padj < p) & (abs(res$log2FoldChange) > log2(fc)))
#   return(res[sig, ])
# }