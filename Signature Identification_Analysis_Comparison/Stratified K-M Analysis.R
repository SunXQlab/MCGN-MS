
###### Install packages in R 
source("http://bioc.ism.ac.jp/biocLite.R")
biocLite("AnnotationDbi")
install.packages("AnnotationDbi")
install.packages("bit")
biocLite("org.Hs.eg.db")
library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)
###### Survival package
install.packages("lattice")
library(lattice)
install.packages("survival")
library(survival)

######## plot KM curves
install.packages("reshape2")
library(reshape2)
install.packages("data.table")
library(data.table)
install.packages("zoo")
library("zoo")
install.packages("survminer")
library("survminer")
install.packages("survival")
library("survival")

#### ROC package
install.packages("gplots")
library(gplots)
install.packages("ROCR")
library(ROCR)
install.packages("dplyr")
library("dplyr")
install.packages("lubridate")
library(lubridate)
install.packages("robustbase")
install.packages("caret")
library(caret)
install.packages("pROC")
library(pROC)

install.packages("timeROC")
install.packages("prodlim")
library(prodlim)
library(quantreg)
install.packages("quantreg")
install.packages("polspline")
library("polspline")
library(timeROC)
install.packages("foreign")
library(foreign)

install.packages("glmnet")
library(glmnet)  # manually load this package by downloading and installing from Toll. 

##################  ensemble ID  to gene symbol   #################################
ensembl2gene <- toTable(org.Hs.egENSEMBL2EG)
gene2symbol <- toTable(org.Hs.egSYMBOL)
ensemble2symbol <- merge(ensembl2gene, gene2symbol, by = 'gene_id')[2:3]
ensemble2symbol=as.matrix(ensemble2symbol)
rownames(ensemble2symbol)=ensemble2symbol[,2]

###TCGA data
setwd("Path")  

Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time')

#####
Gene_GBM=read.csv("TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM,Survival_LGG)


### CGGA data
setwd("Path")   
Gene_CGGA=read.csv("Gene expression_CGGA.csv")
Gene_CGGA=as.matrix(Gene_CGGA)
rownames(Gene_CGGA)=Gene_CGGA[,1]
Gene_CGGA=Gene_CGGA[,-1]

Survival_CGGA=read.csv("patient_survival data_CGGA.csv")
Survival_CGGA=Survival_CGGA[,c(1,3,6,8,9)]
# Survival2=Survival2[Survival2[,3]=='GBM'|Survival2[,3]=='sGBM'|Survival2[,3]=='rGBM',]
rownames(Survival_CGGA)=Survival_CGGA[,1]
dim(Survival_CGGA)
head(Survival_CGGA)



######################## Select genes for signature   #####################################
Marker_genes=c("PXDN" ,  "CDH6", "SCN3A", "ANPEP", "SEMA6B", "CCDC37", "DPP4",  "PRRG1", "GPNMB", "FANCA", "TMEM26", "NETO2" )
# Marker_genes=MNB[A!=0]

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))


Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

PatientNo_1=colnames(Gene_marker)

A=c(0.0111123294, 0.0008619239, -0.8772969022, 0.0016958263, -0.0423078649, 0.0196739563, 0.0013511639, 0.8284922208, 0.0026937362, 1.1843625406, 0.5722500646, 0.1010323338)
# A=A[A!=0]
G=Gene_marker


########  stratified cohert K-M analysis
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

############# Age-stratified coherts #############
Age_LGG=Survival_LGG[,c(1,20)]
rownames(Age_LGG)=Age_LGG[,1]

Age_GBM=Survival_GBM[,c(1,56)]
rownames(Age_GBM)=Age_GBM[,1]                     


Age_TCGA=c(Age_GBM[,2],Age_LGG[,2])
names(Age_TCGA)=c(rownames(Age_GBM),rownames(Age_LGG))
Age=Age_TCGA[PatientNo_1]

##Cohert1
Survival_S=Survival[names(Age[Age<=60]),]
Gene_marker=Gene_marker[,names(Age[Age<=60])]
## Cohert2
Survival_S=Survival[names(Age[Age>60]),]
Gene_marker=Gene_marker[,names(Age[Age>60])]

############# Grade-stratified coherts #############

Grade_TCGA=c(rep('GBM',length(rownames(Gender_GBM))),rep('LGG',length(rownames(Gender_LGG))))
names(Grade_TCGA)=c(rownames(Age_GBM),rownames(Age_LGG))
Grade=Grade_TCGA[PatientNo_1]

Grade[Grade=='LGG']=1
Grade[Grade=='GBM']=2

Grade=as.matrix(t(Grade))

##Cohert1
Survival_S=Survival[names(Grade[,Grade==1]),]
Gene_marker=Gene_marker[,names(Grade[,Grade==1])]

## Cohert2
Survival_S=Survival[names(Grade[,Grade==2]),]
Gene_marker=Gene_marker[,names(Grade[,Grade==2])]

# ################# LDH1 Wild-type/mutation subtype of Glioma patients
# 
# setwd("Path")  
# 
# Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
# rownames(Gene_LGG)=Gene_LGG[,1]
# colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
# Gene_LGG=Gene_LGG[,-1]
# 
# Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
# Survival_LGG=Survival_LGG[,c(1,2,4,67)]
# rownames(Survival_LGG)=Survival_LGG[,1]
# colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','LDH1')
# 
# Survival_LGG=as.matrix(Survival_LGG)
# Survival_LDH1wd=Survival_LGG[Survival_LGG[,4]=='NO',]
# Survival_LDH1mt=Survival_LGG[Survival_LGG[,4]=='YES',]
# 
# Gene_LGG=as.matrix(Gene_LGG)
# Gene_LDH1wd=Gene_LGG[Survival_LGG[,4]=='NO',]
# Gene_LDH1mt=Gene_LGG[Survival_LGG[,4]=='YES',]
# 
# Survival_Gene=Gene_LDH1wd[(Marker_genes),]  # Only 2 gene names were matched, so ...
# colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
# 
# 
# Survival=Survival_LDH1wd[intersect(rownames(Survival_LDH1wd),colnames(Survival_Gene)),]
# Survival=as.matrix(Survival)
# Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]



################ Chemotherapy-stratified cohert ###########################
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,16)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Pharmaceutical_Therapy')

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,62)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','ChemoTherapy')

colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)

Survival_Chemo=rbind(as.matrix(Survival_LGG),as.matrix(Survival_GBM))

##Cohert1
Survival_S=Survival_Chemo[intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='YES',]),rownames(Survival_GBM[Survival_GBM[,4]=='YES',]))),]
Gene_marker=G[,intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='YES',]),rownames(Survival_GBM[Survival_GBM[,4]=='YES',])))]  # G=Gene_marker 
## Cohert2
Survival_S=Survival_Chemo[intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='NO',]),rownames(Survival_GBM[Survival_GBM[,4]=='NO',]))),]
Gene_marker=G[,intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='NO',]),rownames(Survival_GBM[Survival_GBM[,4]=='NO',])))]  # G=Gene_marker 

################ Radiotherapy-stratified cohert ###########################
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,86)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Radio_Therapy')

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,105)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','ChemoTherapy')

colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)

Survival_Radio=rbind(as.matrix(Survival_LGG),as.matrix(Survival_GBM))

##Cohert1
Survival_S=Survival_Chemo[intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='YES',]),rownames(Survival_GBM[Survival_GBM[,4]=='YES',]))),]
Gene_marker=G[,intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='YES',]),rownames(Survival_GBM[Survival_GBM[,4]=='YES',])))]  # G=Gene_marker 
## Cohert2
Survival_S=Survival_Chemo[intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='NO',]),rownames(Survival_GBM[Survival_GBM[,4]=='NO',]))),]
Gene_marker=G[,intersect(colnames(G),c(rownames(Survival_LGG[Survival_LGG[,4]=='NO',]),rownames(Survival_GBM[Survival_GBM[,4]=='NO',])))]  # G=Gene_marker 



################ Stratified Data for K-M analysis #################
status=as.numeric(Survival_S[,2])
time=as.numeric(Survival_S[,3])


time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

# status[status==1&time>365*3]=0   for GBM 

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
status0=as.matrix(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="blue") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((status0),predicted0)
roc0
## optimal combination
opt <- which.max(rowSums(cbind(roc0$sensitivities,roc0$specificities)))
## optimal cut-off point 
sort(predicted0,F)[opt]

##########  K-M survival curves for MNB in TCGA and CGGA datasets
groups=matrix(0,1,length(status))
groups[predicted<=sort(predicted,F)[opt]]=1
groups[predicted>sort(predicted,F)[opt]]=2
# groups[predicted<=median(predicted)]=1
# groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)

# d <- dist(predicted, method = "euclidean") # distance matrix
# fit <- hclust(d, method="ward")
# plot(fit) # display dendogram
# groups <- cutree(fit, k=2) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=2, border="red") 


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))


# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, data = as.data.frame(Gene_marker), linetype = "strata", pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))

# ########### time dependent ROC 
# ROC.DSST<-timeROC(T=time,delta=status,
#                   marker=as.numeric(predicted),cause=1,
#                   weighting="cox",
#                   times=c(3,5,7)*365,ROC=TRUE)
# ROC.DSST
# dev.new()
# plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
# lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
# lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
# lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
# text.legend=c(paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[2],3))),paste("AUC at 7 years:",as.character(round(ROC.DSST$AUC[3],3))))
# legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
# title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
# legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)



##################################################################################################################
##################################################################################################################

### Test datasets
time1=Survival_CGGA[,4]
status1=Survival_CGGA[,5]
status1[is.na(status1)]=0
Genes=rownames(na.omit(Survival_Gene))
# Genes='ORC1'
Gene_marker1= matrix(as.numeric(as.matrix(Gene_CGGA[Genes,])),length(Genes),length(status1))  #

##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status1)]
status0=na.omit(status1)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(status1,as.numeric(predicted1))
AUC=auc(roc1)
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CGGA datasets
groups=matrix(0,1,length(status1))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(time1/365, status1) ~ groups, data = as.data.frame(Gene_marker1))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))


# change line size --> 1
# Change line types by groups (i.e. "strata")
# and change color palette
ggsurvplot(fit,  size = 1,  # change line size
           linetype = "strata", # change line type by groups
           break.time.by = 250, # break time axis by 250
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE # Add p-value
)


# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.y.text.col = TRUE)

####### Time dependent ROC
ROC.DSST<-timeROC(T=time1,delta=status1,
                  marker=as.numeric(predicted1),cause=1,
                  weighting="cox",
                  times=c(1,3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 1 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[2],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[3],3))))
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)


######### Assessment of frequent genetic and genomic alterations in gliomas 

### Test datasets
Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.double(as.matrix(Survival_Gene)),dim(Survival_Gene))
colnames(Gene_marker1)=colnames(Survival_Gene)


time1=time[!is.na(status)]
status1=status[!is.na(status)]
Genes=rownames(na.omit(Survival_Gene))

Gene_marker1=Gene_marker1[,!is.na(status)]
A=c(0.0111123294, 0.0008619239, -0.8772969022, 0.0016958263, -0.0423078649, 0.0196739563, 0.0013511639, 0.8284922208, 0.0026937362, 1.1843625406, 0.5722500646, 0.1010323338)

##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

predicted0=as.numeric(predicted1)

##### IDH1 mutation 
setwd("Path")

Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,67)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','LDH1')

Survival_LGG=as.matrix(Survival_LGG)
Survival_LDH1wd=Survival_LGG[Survival_LGG[,4]=='NO',]
Survival_LDH1mt=Survival_LGG[Survival_LGG[,4]=='YES',]

PatientID=colnames(Gene_marker1)
S1=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_LDH1wd),PatientID)))]
S2=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_LDH1mt),PatientID)))]

Wtest=wilcox.test(S1,S2)  # alternative = c("two.sided", "less", "greater")  "two.sided" (default),
pvalue=Wtest$p.value
pvalue

dataset1 <- data.frame(value = c(S1,S2), group = factor(rep(c("WT", "Mutant"), times = c(length(S1), length(S2)))))

dev.new()
boxplot( value ~ group, dataset1,border = c( "#009E73",'red'),cex = 1,cex.axis=2)  #,col.axis = "#009E73"
 

##### CIMP mythelation 
setwd("Path")

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,27)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','CIMP')

Survival_GBM=as.matrix(Survival_GBM)
Survival_CIMPwd=Survival_GBM[Survival_GBM[,4]=='NON G.CIMP',]
Survival_CIMPmt=Survival_GBM[Survival_GBM[,4]=='G.CIMP',]

PatientID=colnames(Gene_marker1)
S3=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_CIMPwd),PatientID)))]
S4=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_CIMPmt),PatientID)))]

Wtest=wilcox.test(S3,S4)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

dataset2 <- data.frame(value = c(S3,S4), group = factor(rep(c("Unmetylated", "Methylated"), times = c(length(S3), length(S4)))))

dev.new()
boxplot( value ~ group,  notch = F,width = c(0.5,0.5), dataset2,border = c( "#009E73",'red'),cex = 1,cex.axis=2)  #,col.axis = "#009E73"


##### 4 subtypes
setwd("Path")

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,28)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','subtype')

Survival_subtype=Survival_GBM
  
Survival_subtype=as.matrix(Survival_subtype)

Survival_C=Survival_GBM[Survival_GBM[,4]=='Classical',]
Survival_M=Survival_GBM[Survival_GBM[,4]=='Mesenchymal',]
Survival_N=Survival_GBM[Survival_GBM[,4]=='Neural',]
Survival_P=Survival_GBM[Survival_GBM[,4]=='Proneural',]

summary(Survival_GBM[,4])
# Classical Mesenchymal      Neural   Proneural 
# 150         166          90         140 

PatientID=colnames(Gene_marker1)
S_C=predicted0[!is.na(match(PatientID,rownames(Survival_C)))]
S_M=predicted0[!is.na(match(PatientID,rownames(Survival_M)))]
S_N=predicted0[!is.na(match(PatientID,rownames(Survival_N)))]
S_P=predicted0[!is.na(match(PatientID,rownames(Survival_P)))]


Wtest=wilcox.test(S_C,S_M)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue


Wtest=wilcox.test(S_C,S_N)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(S_C,S_P)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(S_M,S_N)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue


Wtest=wilcox.test(S_M,S_P)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

Wtest=wilcox.test(S_N,S_P)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue



dataset <- data.frame(value = c(S_C,S_M,S_N,S_P), group = factor(rep(c("Classical", "Mesenchymal","Neural", "Proneural"), times = c(length(S_C), length(S_M),length(S_N), length(S_P)))))

dev.new()
boxplot( value ~ group,  notch = F, dataset, border = c( "#009E73",'red', 'green','pink'),cex = 1,cex.axis=1.3)  #,col.axis = "#009E73"


######### CGGA set for EGFR and PTEN
##### EGFR mutations 
setwd("Path")

Survival_CGGA=read.csv("patient_survival data_CGGA.csv")
Survival_CGGA=Survival_CGGA[,c(1,8,9,20)]
rownames(Survival_CGGA)=Survival_CGGA[,1]
colnames(Survival_CGGA)=c('SampleID','OS_event','OS_time','EGFR')

Survival_CGGA=as.matrix(Survival_CGGA)
Survival_EGFRwd=Survival_CGGA[Survival_CGGA[,4]==0,]
Survival_EGFRmt=Survival_CGGA[Survival_CGGA[,4]==1,]


MNB=c("PXDN" ,  "CDH6", "SCN3A", "ANPEP", "SEMA6B", "CCDC37", "DPP4",  "PRRG1", "GPNMB", "FANCA", "TMEM26", "NETO2" )
Survival_Gene=Gene_CGGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,3])
time=as.numeric(Survival[,2])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene)),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]


predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted)



PatientID=colnames(Gene_marker)
S5=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_EGFRwd),PatientID)))]
S6=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_EGFRmt),PatientID)))]

Wtest=wilcox.test(S5,S6)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

dataset3 <- data.frame(value = c(S5,S6), group = factor(rep(c("WT", "Mutant"), times = c(length(S5), length(S6)))))

dev.new()
boxplot( value ~ group,  notch = F,width = c(0.5,0.5), dataset3,border = c( "#009E73",'red'),cex = 1,cex.axis=2)  #,col.axis = "#009E73"



##### PTEN mutations 
setwd("Path")

Survival_CGGA=read.csv("patient_survival data_CGGA.csv")
Survival_CGGA=Survival_CGGA[,c(1,8,9,16)]
rownames(Survival_CGGA)=Survival_CGGA[,1]
colnames(Survival_CGGA)=c('SampleID','OS_event','OS_time','PTEN')

Survival_CGGA=as.matrix(Survival_CGGA)
Survival_PTENwd=Survival_CGGA[Survival_CGGA[,4]==0,]
Survival_PTENmt=Survival_CGGA[Survival_CGGA[,4]>0,]


# MNB=c("PXDN" ,  "CDH6", "SCN3A", "ANPEP", "SEMA6B", "CCDC37", "DPP4",  "PRRG1", "GPNMB", "FANCA", "TMEM26", "NETO2" )
# Survival_Gene=Gene_CGGA[t(MNB),]
# colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
# 
# 
# Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
# Survival=as.matrix(Survival)
# Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
# 
# status=as.numeric(Survival[,5])
# time=as.numeric(Survival[,4])
# Genes=rownames(na.omit(Survival_Gene))
# Gene_marker=matrix(as.double(as.matrix(Survival_Gene)),dim(Survival_Gene))
# colnames(Gene_marker)=colnames(Survival_Gene)
# 
# time=time[!is.na(status)]
# Gene_marker=Gene_marker[,!is.na(status)]
# status=status[!is.na(status)]
# 
# 
# predicted=matrix(0,1,dim(Gene_marker)[2])
# for (i in 1:dim(Gene_marker)[2])
# {
#   predicted[i]=sum(A*Gene_marker[,i])
# }
# 
# predicted0=as.numeric(predicted)


PatientID=colnames(Gene_marker)
S7=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_PTENwd),PatientID)))]
S8=predicted0[!is.na(match(PatientID,intersect(rownames(Survival_PTENmt),PatientID)))]

Wtest=wilcox.test(S7,S8)  #"two.sided" (default)
pvalue=Wtest$p.value
pvalue

dataset3 <- data.frame(value = c(S7,S8), group = factor(rep(c("WT", "Mutant"), times = c(length(S7), length(S8)))))

dev.new()
boxplot( value ~ group,  notch = F,width = c(0.5,0.5), dataset3,border = c( "#009E73",'red'),cex = 1,cex.axis=2)  #,col.axis = "#009E73"





