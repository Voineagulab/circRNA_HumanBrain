
setwd("/Volumes/Data1/PROJECTS/circRNAs_final/circRNA_Landscape/DHG_DATA/CompleteDataset/REVISION/")
rm(list=ls())
source("SCRIPTS/Final/Supplementary_Code/Functions.R")
library(splines)
############################################ LOAD DATA ############################################

load("DATA_TABLES/data_DS1.rda");load("DATA_TABLES/data_DS2.rda")
load("DATA_TABLES/ratios_ds1.rda");load("DATA_TABLES/ratios_ds2.rda")
samples1=data_DS1$sampleInfo_DS1; libsizes1=data_DS1$libsizes_DS1
samples2=data_DS2$sampleInfo_DS2; libsizes2=data_DS2$libsizes_DS2
infocolsG=c(1:2); infocolsC=c(1:15); infocolsS=c(1:15)


############################################ Linear Models for DE  ############################################
############################################ Format Covariates
factorcols=c("BroadRegionID",  "Sex",  "SeqBatch", "BrainBank", "Seizures")
samples1$Seizures[which(is.na(samples1$Seizures)==TRUE)]="No"
samples2$Seizures[which(is.na(samples2$Seizures)==TRUE)]="No"
#if (ageAsfactors==TRUE) sample_info$Age=cut(sample_info$Age, breaks=seq(from=0, to=70, by=10))
for (j in which(colnames(samples1)%in%factorcols)) samples1[, j]=as.factor(samples1[,j])
for (j in which(colnames(samples2)%in%factorcols)) samples2[, j]=as.factor(samples2[,j])

####Filtering: circRNAs expressed in at least half of the samples
####DS1

usecirc_ds1=which(rowSums(data_DS1$circData_DS1$circCpm_filter[ , -infocolsC] > 0) > nrow(samples1)/2 )
useData_ds1=list(gene=data_DS1$geneData_DS1$rpkm_filter[,],
                 circ=data_DS1$circData_DS1$circCpm_filter[usecirc_ds1,],
                 sj=data_DS1$sjData_DS1$sjmaxCpm_filter[usecirc_ds1,],
                 ci=ratios_ds1$ci[intersect(which(rowSums(is.na(ratios_ds1$ci)==TRUE) ==0) , usecirc_ds1) ,])

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
design.complete_ds1="Neurons + Sex + SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL" # all cov (neuronal estimates based on Zhang et al.)
design.complete_ds1_F5="F5_Neurons + Sex + SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL"# all cov (neuronal estimates based on FANTOM5)
design.reduced_ds1="Sex + SeqBatch + BrainBank + age.splines_ds1 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL" # all cov except cell comp
design.complete_ds2="Neurons + Sex + age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL" # all cov (neuronal estimates based on Zhang et al.)
design.complete_ds2_F5="F5_Neurons + Sex + age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL" # all cov (neuronal estimates based on FANTOM5)
design.reduced_ds2="Sex + age.splines_ds2 + BroadRegionID + Median.3prime.Bias.picard + Seizures + ASD.CTL" # all cov except cell comp
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

dat2=useData_ds2[[j]]
model.complete_ds2[[j]]=lm(paste0("dat2 ~", design.complete_ds2),  samples2)
model.complete_ds2_F5[[j]]=lm(paste0("dat2 ~", design.complete_ds2_F5),  samples2)
model.reduced_ds2[[j]]=lm(paste0("dat2 ~", design.reduced_ds2),  samples2)

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

varnames=c("age.splines_1","age.splines_2", "age.splines_3", "age.splines_4","BroadRegionIDCTX", "SeizuresYes", "Median.3prime.Bias.picard", "ASD.CTLCTL")
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

DEresults_AllSamples=list(results=results, results.summary=results.summary)
save(DEresults_AllSamples,file="Analysis/DE/DEresults_AllSamples.rda")

