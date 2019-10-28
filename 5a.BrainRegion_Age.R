
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
rm(list=ls())
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
library(splines)
par(pch=20, col="grey", family="serif")
############################################ LOAD DATA ############################################

load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/ratios_ds1.rda");load("DATA_TABLES/ratios_ds2.rda")
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)

############################################ PCA Plots ############################################
pdf("Analysis/DatasetCharacterisationPlots/Fig5.PCA_Plots.pdf", height=4, width=8)
par(mfrow=c(1,2))
plotPCA(ratios_ds1$clr , infocolsC, samples1, "CLR,DS1", obs = "complete")
plotPCA(ratios_ds2$clr  , infocolsC, samples2, "CLR, DS2", obs = "complete")
dev.off()

pdf("Analysis/DatasetCharacterisationPlots/SuppFig.PCAplots_all.pdf", height=8, width=8)
par(mfrow=c(3,3))
plotPCA(data_DS1$geneData_DS1$rpkm_filter, infocolsG, samples1, "Gene Expression,  DS1", obs = "complete")
plotPCA(data_DS1$circData_DS1$circCpm_filter , infocolsC, samples1, "CircRNAs,  DS1", obs = "complete")
plotPCA(ratios_ds1$ci , infocolsC, samples1, "CI,  DS1", obs = "complete")

plotPCA(data_DS2$geneData_DS2$rpkm_filter , infocolsG, samples2, "Gene Expression,DS2", obs = "complete")
plotPCA(data_DS2$circData_DS2$circCpm_filter , infocolsC, samples2, "CircRNAs,  DS2", obs = "complete")
plotPCA(ratios_ds2$ci , infocolsC, samples2, "CI,  DS2", obs = "complete")
dev.off()
#revision: plot perc of variance explained
# data= ratios_ds1$clr ; infocols= infocolsC; samples= samples1; title= "CLR,DS1"; obs = "complete"
# samples=samples[match(colnames(data)[-infocols], samples$Sample) , ]
# colors <- factor(samples$RegionID, levels=(c("ba9", "ba41-42-22", "vermis")) ,labels=c("#66c2a5",  "#8da0cb" ,"#fc8d62"))
# p=rep(17, nrow(samples))
# p[which(samples$ASD.CTL %in%"ASD")]=16
# cormx=cor(data[,-infocols], method="s", use=obs)
# pca=princomp(cormx)
# plot(pca$loadings[,1], pca$loadings[,2], pch=p, col=as.character(colors),  main=title)
# var= pca$sdev^2
# percvar=var/ sum(var)
# barplot(percvar[1:10])


############################################ Linear Models for DE  ############################################
############################################ Format Covariates
factorcols=c("BroadRegionID",  "Sex",  "SeqBatch", "BrainBank")
#if (ageAsfactors==TRUE) sample_info$Age=cut(sample_info$Age, breaks=seq(from=0, to=70, by=10))
for (j in which(colnames(samples1)%in%factorcols)) samples1[, j]=as.factor(samples1[,j])
for (j in which(colnames(samples2)%in%factorcols)) samples2[, j]=as.factor(samples2[,j])

####Filtering: circRNAs expressed in at least half of the CTL samples
####DS1
con_ds1=as.character(samples1$Sample[which(samples1$ASD.CTL == "CTL")])
usecirc_ds1=which(rowSums(data_DS1$circData_DS1$circCpm_filter[ , con_ds1 ] > 0) > length(con_ds1)/2 )
useData_ds1=list(gene=data_DS1$geneData_DS1$rpkm_filter[,con_ds1],
                 circ=data_DS1$circData_DS1$circCpm_filter[usecirc_ds1,con_ds1],
                 sj=data_DS1$sjData_DS1$sjmaxCpm_filter[usecirc_ds1,con_ds1],
                 ci=ratios_ds1$ci[intersect(which(rowSums(is.na(ratios_ds1$ci)==TRUE) ==0) , usecirc_ds1) ,con_ds1])
libsizes1=libsizes1[con_ds1];samples1=samples1[con_ds1,]
####DS2
usecirc_ds2=which(rowSums(data_DS2$circData_DS2$circCpm_filter[ , -infocolsC] > 0) > nrow(samples2)/2 )
useData_ds2=list(gene=data_DS2$geneData_DS2$rpkm_filter,
                   circ=data_DS2$circData_DS2$circCpm_filter[usecirc_ds2,],
                   sj=data_DS2$sjData_DS2$sjmaxCpm_filter[usecirc_ds2,],
                   ci=ratios_ds2$ci[intersect(which(rowSums(is.na(ratios_ds2$ci)==TRUE) ==0) , usecirc_ds2) ,])

age.splines_ds1=bs(samples1$Age, knots=c(0,10,20,40), degree=1)
age.splines_ds2=bs(samples2$Age, knots=c(0,10,20,40), degree=1)
##############Designs
#Note: same analysis was done with Age as numeric variable
design.complete_ds1="Neurons + Sex +  SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard" # all cov (neuronal estimates based on Zhang et al.)
design.complete_ds1_F5="F5_Neurons + Sex + SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard "# all cov (neuronal estimates based on FANTOM5)
design.reduced_ds1="Sex + RIN + SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard" # all cov except cell comp
design.complete_ds2="Neurons + Sex + age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard" # all cov (neuronal estimates based on Zhang et al.)
design.complete_ds2_F5="F5_Neurons + Sex + age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard" # all cov (neuronal estimates based on FANTOM5)
design.reduced_ds2="Sex +  age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard" # all cov except cell comp
##############Format data
useData_ds1=lapply(useData_ds1, format_for_lm, libs=libsizes1, samples=samples1, voom=FALSE)
useData_ds2=lapply(useData_ds2, format_for_lm, libs=libsizes2, samples=samples2, voom=FALSE)
##############Models
model.complete_ds1=list();model.complete_ds1_F5=list();model.reduced_ds1=list()
model.complete_ds2=list();model.complete_ds2_F5=list();model.reduced_ds2=list()
for (  j in c(1:4))
{
dat1=useData_ds1[[j]]
model.complete_ds1[[j]]=lm(paste0("dat1 ~", design.complete_ds1),  samples1)
model.complete_ds1_F5[[j]]=lm(paste0("dat1 ~", design.complete_ds1_F5),  samples1)
model.reduced_ds1[[j]]=lm(paste0("dat1 ~", design.reduced_ds1),  samples1)
print (paste("Done", names(useData_ds1)[j]))
dat2=useData_ds2[[j]]
model.complete_ds2[[j]]=lm(paste0("dat2 ~", design.complete_ds2),  samples2)
model.complete_ds2_F5[[j]]=lm(paste0("dat2 ~", design.complete_ds2_F5),  samples2)
model.reduced_ds2[[j]]=lm(paste0("dat2 ~", design.reduced_ds2),  samples2)
print (paste("Done", names(useData_ds2)[j]))
names(model.complete_ds1)[[j]]=names(model.complete_ds1_F5)[[j]]=names(model.reduced_ds1)[[j]]=names(useData_ds1)[[j]]
names(model.complete_ds2)[[j]]=names(model.complete_ds2_F5)[[j]]=names(model.reduced_ds2)[[j]]=names(useData_ds2)[[j]]

}

##############Pvals with multiple testing adjustment
pvals=list()
pvals$model.complete_ds1=lapply(model.complete_ds1, 
                            function( dat) adjPvals.from.lm(coefficients(summary(dat))))
pvals$model.complete_ds2=lapply(model.complete_ds2, 
                                function( dat) adjPvals.from.lm(coefficients(summary(dat))))
pvals$model.complete_ds1_F5=lapply(model.complete_ds1_F5, 
                                function( dat) adjPvals.from.lm(coefficients(summary(dat))))
pvals$model.complete_ds2_F5=lapply(model.complete_ds2_F5, 
                                function( dat) adjPvals.from.lm(coefficients(summary(dat))))
pvals$model.reduced_ds1=lapply(model.reduced_ds1, 
                                function( dat) adjPvals.from.lm(coefficients(summary(dat))))
pvals$model.reduced_ds2=lapply(model.reduced_ds2, 
                               function( dat) adjPvals.from.lm(coefficients(summary(dat))))

for (j in c(1:6)) pvals[[j]]=lapply(pvals[[j]], 
                                    function(x){
                                      colnames(x)=gsub("Pr(>|t|)_", "", colnames(x), fixed=TRUE);
                                      colnames(x)=gsub("(", "", colnames(x), fixed=TRUE);
                                      colnames(x)=gsub(")", "", colnames(x), fixed=TRUE);
                                      colnames(x)=gsub("ds1", "", colnames(x), fixed=TRUE);
                                      colnames(x)=gsub("ds2", "", colnames(x), fixed=TRUE);
                                      return(x)
                                    })
varnames=c("age.splines_1","age.splines_2", "age.splines_3", "age.splines_4","BroadRegionIDCTX")
for (j in c(1:6)) pvals[[j]]=lapply(pvals[[j]], 
                                    function(x) x=x[, which(colnames(x)%in%varnames)] )

##############Residuals
residuals=list()
residuals$model.complete_ds1=lapply(model.complete_ds1, 
                                    function(dat) t(residuals(dat)))
residuals$model.complete_ds2=lapply(model.complete_ds2, 
                                    function(dat) t(residuals(dat)))
residuals$model.complete_ds1_F5=lapply(model.complete_ds1_F5, 
                                       function(dat) t(residuals(dat)))
residuals$model.complete_ds2_F5=lapply(model.complete_ds2_F5, 
                                       function(dat) t(residuals(dat)))
residuals$model.reduced_ds1=lapply(model.reduced_ds1, 
                                   function(dat) t(residuals(dat)))
residuals$model.reduced_ds2=lapply(model.reduced_ds2, 
                                   function(dat) t(residuals(dat)))

##############Significant changes
results=list()
results.summary=list()
#######Looping through all contrasts
for(k in c(1:length(varnames)))
{#######Looping through all models
  print("=====");print(varnames[k])
  var=list()
  for (j in c(1:length(pvals))) 
  {var[[j]]= lapply(pvals[[j]], function(dat) rownames(dat)[which(dat[,k] < 0.05 )] )
  print(j)
  }
  names(var)=names(pvals)
  results[[k]]=var
  names(results)[k]=varnames[k]
  
  var.summary=as.data.frame(matrix(NA, nrow=9, ncol=4)); 
  colnames(var.summary)=names(var[[j]])
  rownames(var.summary)=c(names(var), "both_ds_model.complete", "both_ds_model.complete_F5","both_ds_model.reduced")
  for( j in c(1:length(var)))
    var.summary[j,]=lapply (var[[j]] , function (x) length(x))   
  for( j in c(1:length(var[[1]])))
  {
    var.summary["both_ds_model.complete", j]=length(intersect(var$model.complete_ds1[[j]], var$model.complete_ds2[[j]]))
    var.summary["both_ds_model.complete_F5", j]=length(intersect(var$model.complete_ds1_F5[[j]], var$model.complete_ds2_F5[[j]]))
    var.summary["both_ds_model.reduced", j]=length(intersect(var$model.reduced_ds1[[j]], var$model.reduced_ds2[[j]]))
  }
  results.summary[[k]]=var.summary
  names(results.summary)[k]=varnames[k]
  print(var.summary)
}
# Save all LM results
DEresults_CTL=list(results=results, results.summary=results.summary)
save(DEresults_CTL,file="Analysis/DE/DEresults_CTL.rda")

# Save and plot data for Region Significant circRNAs
ci.ds1=ratios_ds1$ci; ci.ds2=ratios_ds2$ci
signif.circ=intersect(results$BroadRegionIDCTX$model.complete_ds1$ci, results$BroadRegionIDCTX$model.complete_ds2$ci)
signif.circ.data1=ci.ds1[signif.circ, ]
signif.circ.data2=ci.ds2[signif.circ, ]
signif.circ.info=signif.circ.data1[,infocolsC]
signif.circ.info$pval.ds1=pvals$model.complete_ds1$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds1$ci))]
signif.circ.info$pval.ds2=pvals$model.complete_ds2$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds2$ci))]
signif.circ.info$pval.ds1_F5=pvals$model.complete_ds1_F5$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds1_F5$ci))]
signif.circ.info$pval.ds2_F5=pvals$model.complete_ds2_F5$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds2_F5$ci))]
signif.circ.info$meanCI.CTX.ds1=apply(signif.circ.data1[, grep("ba", colnames(signif.circ.data1))], 1, mean)
signif.circ.info$meanCI.CB.ds1=apply(signif.circ.data1[, grep("vermis", colnames(signif.circ.data1))], 1, mean)
signif.circ.info$meanCI.CTX.ds2=apply(signif.circ.data2[, grep("ba", colnames(signif.circ.data2))], 1, mean)
signif.circ.info$meanCI.CB.ds2=apply(signif.circ.data2[, grep("vermis", colnames(signif.circ.data2))], 1, mean)
signif.circ.info$meanCI.Diff.ds1=signif.circ.info$meanCI.CB.ds1-signif.circ.info$meanCI.CTX.ds1
signif.circ.info$meanCI.Diff.ds2=signif.circ.info$meanCI.CB.ds2-signif.circ.info$meanCI.CTX.ds2
write.csv(signif.circ.info, "Analysis/DE/SuppTable.region_circ.csv")
pdf("Analysis/DE/Region.meanCI.Diff.pdf", height=6, width=6)
s=order(abs(signif.circ.info$meanCI.Diff.ds1), decreasing = TRUE)
p=c(which(signif.circ.info$meanCI.Diff.ds1 < 0), s[1:20])
plot(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2,
     sub=round(cor(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2, method="s"),2),
     pch=20,
     col=rgb(0,0,0,0.5),
     xlab="Mean CI difference in DS1",
     ylab="Mean CI difference in DS2")
text(signif.circ.info$meanCI.Diff.ds1[p], 
     signif.circ.info$meanCI.Diff.ds2[p],
     signif.circ.info$Symbol[p],
     cex=0.3,pos=3)
dev.off()
pdf("Analysis/DE/Region.meanCI.Diff.large.pdf", height=15, width=15)
s=order(abs(signif.circ.info$meanCI.Diff.ds1), decreasing = TRUE)
p=c(which(signif.circ.info$meanCI.Diff.ds1 < 0), s[1:20])
plot(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2,
     sub=round(cor(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2, method="s"),2),
     pch=20,
     col=rgb(0,0,0,0.5),
     xlab="Mean CI difference in DS1",
     ylab="Mean CI difference in DS2")
text(signif.circ.info$meanCI.Diff.ds1[p], 
     signif.circ.info$meanCI.Diff.ds2[p],
     signif.circ.info$Symbol[p],
     cex=0.3,pos=3)
dev.off()

length(results$BroadRegionIDCTX$model.complete_ds1$ci)
# 496
length(intersect(results$BroadRegionIDCTX$model.complete_ds1$ci, results$BroadRegionIDCTX$model.complete_ds2$ci))
#[1] 205
length(intersect(results$BroadRegionIDCTX$model.reduced_ds1$ci, results$BroadRegionIDCTX$model.reduced_ds2$ci))
#[1] 408
length(intersect(intersect(results$BroadRegionIDCTX$model.reduced_ds1$ci, results$BroadRegionIDCTX$model.reduced_ds2$ci),intersect(results$BroadRegionIDCTX$model.complete_ds1$ci, results$BroadRegionIDCTX$model.complete_ds2$ci)))
#[1] 197

##Save and plot data for Region Significant circRNAs - Without correction for proportion of neurons
signif.circ=intersect(results$BroadRegionIDCTX$model.reduced_ds1$ci, results$BroadRegionIDCTX$model.reduced_ds2$ci)
signif.circ.data1=ci.ds1[signif.circ, ]
signif.circ.data2=ci.ds2[signif.circ, ]
signif.circ.info=signif.circ.data1[,infocolsC]
signif.circ.info$pval.ds1=pvals$model.complete_ds1$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds1$ci))]
signif.circ.info$pval.ds2=pvals$model.complete_ds2$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds2$ci))]
signif.circ.info$pval.ds1_F5=pvals$model.complete_ds1_F5$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds1_F5$ci))]
signif.circ.info$pval.ds2_F5=pvals$model.complete_ds2_F5$ci$BroadRegionIDCTX[match(signif.circ, rownames(pvals$model.complete_ds2_F5$ci))]
signif.circ.info$meanCI.CTX.ds1=apply(signif.circ.data1[, grep("ba", colnames(signif.circ.data1))], 1, mean)
signif.circ.info$meanCI.CB.ds1=apply(signif.circ.data1[, grep("vermis", colnames(signif.circ.data1))], 1, mean)
signif.circ.info$meanCI.CTX.ds2=apply(signif.circ.data2[, grep("ba", colnames(signif.circ.data2))], 1, mean)
signif.circ.info$meanCI.CB.ds2=apply(signif.circ.data2[, grep("vermis", colnames(signif.circ.data2))], 1, mean)
signif.circ.info$meanCI.Diff.ds1=signif.circ.info$meanCI.CB.ds1-signif.circ.info$meanCI.CTX.ds1
signif.circ.info$meanCI.Diff.ds2=signif.circ.info$meanCI.CB.ds2-signif.circ.info$meanCI.CTX.ds2
write.csv(signif.circ.info, "Analysis/DE/SuppTable.region_circ_noCellTypeCorrection.csv")
pdf("Analysis/DE/Region.meanCI.Diff_noCellTypeCorrection.pdf", height=3, width=5)
plot(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2,sub=cor(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2, method="s"))
dev.off()
pdf("Analysis/DE/Region.meanCI.Diff.large_noCellTypeCorrection.pdf", height=10, width=20)
plot(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2,sub=cor(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2, method="s"))
text(signif.circ.info$meanCI.Diff.ds1, signif.circ.info$meanCI.Diff.ds2,signif.circ.info$Symbol )
dev.off()


# nc=setdiff(intersect(results$BroadRegionIDCTX$model.reduced_ds1$ci, results$BroadRegionIDCTX$model.reduced_ds2$ci),intersect(results$BroadRegionIDCTX$model.complete_ds1$ci, results$BroadRegionIDCTX$model.complete_ds2$ci))
# cor_nc_1=cor(t(signif.circ.data1[nc,rownames(samples1)]), samples1$Neurons, method="s")
# plot(cor_nc_1,signif.circ.info[nc, "meanCI.Diff.ds1"])
# 
# cor_nc_1=cor(t(signif.circ.data1[nc,rownames(samples1)]), samples1$Neurons, method="s")
# plot(cor_nc_1,signif.circ.info[nc, "meanCI.Diff.ds1"])
############################################ N circ RNAs expressed vs Age
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)

circ_ds1=data_DS1$circData_DS1$circCpm_filter
circ_ds2=data_DS2$circData_DS2$circCpm_filter
samples1$nCirc1=colSums(circ_ds1[, -infocolsC] > 0 , na.rm=TRUE)
samples2$nCirc2=colSums(circ_ds2[, -infocolsC] > 0 , na.rm=TRUE)

ctx=which(samples1$BroadRegionID%in%"CTX")
cb=which(samples1$BroadRegionID%in%"CB")
CTL=which(samples1$ASD.CTL%in%"CTL")

pdf("Analysis/DE/Fig5.Age.Ncirc.pdf", height=3, width=5)
plotC <- ggplot(samples1[CTL,], aes(y = nCirc1, x = Age, col = BroadRegionID)) +
  geom_point(aes(shape = BroadRegionID), size = 1.5) +
  geom_smooth(fill="lightgrey") +
  scale_colour_manual(values=c("#fc8d62","#66c2a5" )) +
  ylab ("Number of circRNAs expressed")+
  theme(axis.text = element_text(family="serif",size=10),
        axis.title = element_text(family="serif",size=10),
        legend.text = element_text(family="serif",size=10))
plotC
dev.off()
summary(lm(nCirc1~ Neurons + Sex + Median.3prime.Bias.picard + SeqBatch + BrainBank +Age, samples1[CTL,]))
summary(lm(nCirc2~ Neurons + Sex + Median.3prime.Bias.picard + Age, samples2[which(samples2$ASD.CTL%in%"CTL"),]))
# Call:
#   lm(formula = nCirc1 ~ Neurons + Sex + Median.3prime.Bias.picard + 
#        SeqBatch + BrainBank + Age, data = samples1[CTL, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2448.9  -793.2    -8.4   670.9  3343.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -15356.46    4696.93  -3.269 0.001787 ** 
#   Neurons                    -2103.87     983.47  -2.139 0.036496 *  
#   SexM                         443.77     372.24   1.192 0.237892    
# Median.3prime.Bias.picard  24138.98    6542.39   3.690 0.000486 ***
#   SeqBatchbatch2             -1517.37     615.99  -2.463 0.016651 *  
#   SeqBatchbatch3             -2525.65     710.12  -3.557 0.000741 ***
#   BrainBankNICHD                30.29     448.03   0.068 0.946330    
# Age                           23.99      12.82   1.871 0.066279 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1182 on 60 degrees of freedom
# Multiple R-squared:  0.3569,	Adjusted R-squared:  0.2818 
# F-statistic: 4.756 on 7 and 60 DF,  p-value: 0.0002669
