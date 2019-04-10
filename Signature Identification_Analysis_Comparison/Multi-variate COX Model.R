
###### Install packages in R 
source("http://bioc.ism.ac.jp/biocLite.R")
biocLite("AnnotationDbi")
install.packages("AnnotationDbi")
install.packages("bit")
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
Gene_CCGA=read.csv("Gene expression_CGGA.csv")
Gene_CCGA=as.matrix(Gene_CCGA)
rownames(Gene_CCGA)=Gene_CCGA[,1]
Gene_CCGA=Gene_CCGA[,-1]

Survival_CCGA=read.csv("patient_survival data_CGGA.csv")
Survival_CCGA=Survival_CCGA[,c(1,3,6,8,9)]
# Survival2=Survival2[Survival2[,3]=='GBM'|Survival2[,3]=='sGBM'|Survival2[,3]=='rGBM',]
rownames(Survival_CCGA)=Survival_CCGA[,1]
dim(Survival_CCGA)
head(Survival_CCGA)



######################## Select genes for signature   #####################################
# DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
# DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
# DNG=t(as.matrix(DN_human[,1]))
# MNB=as.character(DNG)  # Macrophage-associated gene signature

##################################################################################################################
##################################################################################################################

### Evaluate Multivariate COX model on Test dataset  (TCGA set)
Survival_Gene=Gene_TCGA[t(MNB),]
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

### COX Model
# fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")

surv=Surv(as.double(time),as.double(status))

# cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1)
# dev.new()
# plot(cvfit)
# A=coef(cvfit, s = "lambda.min")
# which(A!=0)
# sum(A!=0)
# MNB[which(A!=0)]
# 
# A0=A

# A=load('Model_Coefficient.RDatac')
########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival status

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A0*Gene_marker[,i])
}

# predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# status0=as.numeric(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="red") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# 
# roc0=roc(status0,predicted0)
# AUC=auc(roc0)
# AUC

troc0<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_1=predicted

##########  EGFR
Marker_genes=c('EGFR')  # R+TFs+Liganfs

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
names(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[!is.na(status)]
status=status[!is.na(status)]

PatientNo_2=names(Gene_marker)
### COX Model
fit <- coxph(Surv(time, status) ~ (Gene_marker), method="breslow")
fit
cox.zph(fit)
A1=fit$coefficients

########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,length(Gene_marker))
for (i in 1:length(Gene_marker))
{
  predicted[i]=sum(A1*Gene_marker[i])
}

predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# status0=as.numeric(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new = TRUE)
# plot(perf,colorize=FALSE, col="blue") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# 
# roc1=roc(status0,predicted0)
# AUC=auc(roc1)
# AUC

troc1<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_2=predicted0



##########  Cheng et al Neurology 2016; 86 (24)  
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112每118


Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

PatientNo_4=colnames(Gene_marker)

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
fit
cox.zph(fit)
A3=fit$coefficients

A3=c(-0.6718,0.1658, 0.2584, -0.1811, 0.1165, -0.4046, 0.1543, 0.1223)  # values in that paper

########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A3*Gene_marker[,i])
}

# predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# status0=as.numeric(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new = TRUE)
# plot(perf,colorize=FALSE, col="yellow") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# 
# roc3=roc(status0,predicted0)
# AUC=auc(roc3)
# AUC

troc3<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_3=predicted



### IGF1-PI3K Pathway (Signature 4) 
################################################################
Marker_genes=c('PIK3R1',	'PIK3R2'	,'PIK3R3',	'PIK3R4',	'PIK3R5'	,'PIK3R6',	'PIK3AP1',		'AKT1',	'AKT2',	'AKT3',		'IGF1',	'IGF1R',		'IL4',	'IL4R',		'NFATC1',	'NFATC2',	'NFATC3',	'NFATC4'	,'NFAT5',		'STAT6','CSF1','CSF1R')
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGA))
Marker_genes=intersect(Marker_genes,rownames(Gene_TCGA))


Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

PatientNo_4=colnames(Gene_marker)


A=A4
# A4
# [1] -0.08057491  0.15705491  0.10197133 -0.03924833  6.38801401  0.27136837  0.14195919  0.01471263  0.39661697 -0.05663284  1.52427987
# [12]  0.02936325 -5.81479105 -0.06172801 -7.74822977 11.42249530 -0.23941925 -5.91291532 -1.91197970  0.04755584 -0.02362737


########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A4*Gene_marker[,i])
}


troc3<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_4=predicted


######### Age, Grade, Gender

## Read clinical data
setwd("Path")  
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")

## Age
Age_LGG=Survival_LGG[,c(1,20)]
rownames(Age_LGG)=Age_LGG[,1]

Age_GBM=Survival_GBM[,c(1,56)]
rownames(Age_GBM)=Age_GBM[,1]                     
 

Age_TCGA=c(Age_GBM[,2],Age_LGG[,2])
names(Age_TCGA)=c(rownames(Age_GBM),rownames(Age_LGG))
Age=Age_TCGA[PatientNo_1]

Age[Age<=60]=1
Age[Age>60]=2
Age=as.matrix(t(Age))

# Gender
Gender_LGG=Survival_LGG[,c(1,50)]
rownames(Gender_LGG)=Gender_LGG[,1]

Gender_GBM=Survival_GBM[,c(1,76)]
rownames(Gender_GBM)=Gender_GBM[,1]  

Gender_TCGA=c(as.vector(Gender_GBM[,2]), as.vector(Gender_LGG[,2]))
names(Gender_TCGA)=c(rownames(Gender_GBM),rownames(Gender_LGG))
Gender=Gender_TCGA[PatientNo_1]

Gender[Gender=='FEMALE']=1
Gender[Gender=='MALE']=2
Gender=as.matrix(t(Gender))


# Grade

Grade_TCGA=c(rep('GBM',length(rownames(Gender_GBM))),rep('LGG',length(rownames(Gender_LGG))))
names(Grade_TCGA)=c(rownames(Gender_GBM),rownames(Gender_LGG))
Grade=Grade_TCGA[PatientNo_1]
  
Grade[Grade=='LGG']=1
Grade[Grade=='GBM']=2

Grade=as.matrix(t(Grade))

################# Multi-variable COX model
Marker=rbind(as.numeric(Age),as.numeric(Gender),as.numeric(Grade),Score_1,as.matrix(t(Score_2)),Score_3,Score_4)

fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")

summary(fit)

################# Uni-variable COX model
Marker=rbind(as.numeric(Age))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Gender))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Grade))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)


Marker=rbind(Score_1)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(Score_2)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 
 
Marker=rbind(Score_3)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)  
 
Marker=rbind(Score_4)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 

######## 5-year Survival evaluation

status[status==1&time>365*5]=0  # 5-year survival status

################# Multi-variable COX model
Marker=rbind(as.numeric(Age),as.numeric(Gender),as.numeric(Grade),Score_1,as.matrix(t(Score_2)),Score_3,Score_4)

fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")

summary(fit)

################# Uni-variable COX model
Marker=rbind(as.numeric(Age))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Gender))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Grade))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)


Marker=rbind(Score_1)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(Score_2)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 

Marker=rbind(Score_3)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)  

Marker=rbind(Score_4)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 






################ LASSO regression for multivariate COX model
Marker=rbind(as.numeric(Age),as.numeric(Gender),as.numeric(Grade),Score_1,as.matrix(t(Score_2)),Score_3)

surv=Surv(as.double(time),as.double(status))

cvfit<-cv.glmnet(t(Marker),surv,family="cox",alpha=1)
plot(cvfit)
AA=as.numeric(coef(cvfit, s = "lambda.min"))
# A=coef(cvfit, s="lambda.1se")
which(AA!=0)
sum(AA!=0)
# Marker[which(AA!=0),]

fit <- coxph(Surv(time, status) ~ t(Marker[which(AA!=0),]), method="breslow")
summary(fit)

# scoxm=step(fit)
# 
# ###### Nomograme 
# 
# install.packages("rms")
# library(rms)
# require(rms)








##################################################################################################################
### Training datasets  (CGGA dataset)
##################################################################################################################



### MNB
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG)  # Macrophage-associated gene signature

# Survival_Gene=Gene_TCGA[t(MNB),]
# colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
# Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
# Survival=as.matrix(Survival)
# Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
# Genes=rownames(na.omit(Survival_Gene))

Genes=MNB

Gene_marker1=matrix(as.numeric(as.matrix(Gene_CCGA[Genes,])),dim(Gene_CCGA[Genes,]))


A=A0

##### ROC for test set
time1=Survival_CCGA[,4]
status1=Survival_CCGA[,5]
status1[is.na(status1)]=0
# status1[status1==1 & time1>365*5]=0  # 5-year survival

predicted1=matrix(0,1,dim(Gene_marker1)[2])
for (i in 1:dim(Gene_marker1)[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

# predicted0=as.numeric(predicted1)
# predicted0=predicted1[!is.na(status1)]
# status0=na.omit(status1)
# pred <- prediction(as.numeric(predicted0),status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="red") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# roc4=roc(status0,as.numeric(predicted1),auc.polygon=TRUE)
# roc4
Score_1=predicted1

###########################
###################### time-dependent ROC  
troc4<-timeROC(T=time1,delta=status1,
               marker=as.numeric(predicted1),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)


###############################
###############################


##########  EGFR
################################################################
Marker_genes=c('EGFR')  # R+TFs+Liganfs

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

time1=Survival_CCGA[,4]
status1=Survival_CCGA[,5]
status1[is.na(status1)]=0
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.numeric(as.matrix(Gene_CCGA[Genes,])),dim(as.matrix(Gene_CCGA[Genes,])))


A=A1
##### ROC for test set
# status1[status1==1&time1>365*5]=0  # 5-year survival
Gene_marker1=as.matrix(t(Gene_marker1))
predicted1=matrix(0,1,dim(Gene_marker1)[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

# predicted0=as.numeric(predicted1)
# predicted0=predicted1[!is.na(status1)]
# status0=na.omit(status1)
# pred <- prediction(as.numeric(predicted0),status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new = TRUE)
# plot(perf,colorize=FALSE, col="blue") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# roc5=roc(status1,as.numeric(predicted1))
# AUC=auc(roc5)
# AUC


troc5<-timeROC(T=time1,delta=status1,
               marker=as.numeric(predicted1),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_2=predicted1

# ### Bao et al.
# ################################################################
# Marker_genes=c('BIRC5', 'TEAD2', 'TUBA1B', 'MT1E','RAB1A','SFXN4','TPX2','HDAC4','FAM125B')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112每118
# ## "FAM125B" is also named as "MVB12B"
# 
# # Marker_genes=rownames(Gene_TCGA)  # select genes from all of TCGA genes
# 
# 
# Survival_Gene=Gene_TCGA[t(Marker_genes),]
# colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
# Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
# Survival=as.matrix(Survival)
# Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
# 
# time1=Survival_CCGA[,4]
# status1=Survival_CCGA[,5]
# status1[is.na(status1)]=0
# Genes=rownames(na.omit(Survival_Gene))
# Genes=intersect(Genes,rownames(Gene_CCGA))
# Genes=c('BIRC5', 'TEAD2', 'TUBA1B', 'MT1E','RAB1A','SFXN4','TPX2','HDAC4','MVB12B')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112每118
# ## "FAM125B" is also named as "MVB12B"
# Gene_marker1=matrix(as.numeric(as.matrix(Gene_CCGA[Genes,])),dim(Gene_CCGA[Genes,]))
# 
# A=A2[!is.na(match(Marker_genes,Genes))]
# 
# A=c(1.653,2.382, 2.942, 1.365, 0.584, 0.22, 2.588, 0.548, 2.447)
# A=log2(A)
# 
# ##### ROC for test set
# # status1[status1==1&time1>365*5]=0  # 5-year survival
# 
# predicted1=matrix(0,1,dim(Gene_marker1)[2])
# for (i in 1:dim(Gene_marker1)[2])
# {
#   predicted1[i]=sum(A*Gene_marker1[,i])
# }
# 
# # predicted0=as.numeric(predicted1)
# # predicted0=predicted1[!is.na(status1)]
# # status0=na.omit(status1)
# # pred <- prediction(as.numeric(predicted0),status0)
# # perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# # performance(pred,"auc") # shows calculated AUC for model
# # par(new = TRUE)
# # plot(perf,colorize=FALSE, col="green") # plot ROC curve
# # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# # 
# # roc6=roc(status1,as.numeric(predicted1))
# # AUC=auc(roc6)
# # AUC
# # 
# troc6<-timeROC(T=time1,delta=status1,
#                marker=as.numeric(predicted1),cause=1,
#                iid=TRUE,
#                times=c(3,5)*365,ROC=TRUE)
# 
# 
# Score_3=predicted1

### Cheng et al. 
################################################################
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112每118

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

time1=Survival_CCGA[,4]
status1=Survival_CCGA[,5]
status1[is.na(status1)]=0
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.numeric(as.matrix(Gene_CCGA[Genes,])),dim(Gene_CCGA[Genes,]))

A=A3
# 
A=c(-0.6718,0.1658, 0.2584, -0.1811, 0.1165, -0.4046, 0.1543, 0.1223)

#### ROC for test set

# status1[status1==1&time1>365*5]=0  # 5-year survival

predicted1=matrix(0,1,dim(Gene_marker1)[2])
for (i in 1:dim(Gene_marker1)[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

# predicted0=as.numeric(predicted1)
# predicted0=predicted1[!is.na(status1)]
# status0=na.omit(status1)
# pred <- prediction(as.numeric(predicted0),status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new = TRUE)
# plot(perf,colorize=FALSE, col="yellow") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )

# roc7=roc(status1,as.numeric(predicted1))
# AUC=auc(roc7)
# AUC


troc7<-timeROC(T=time1,delta=status1,
               marker=as.numeric(predicted1),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_3=predicted1


### IGF1-PI3K Pathway (Signature 4) 
################################################################
Marker_genes=c('PIK3R1',	'PIK3R2'	,'PIK3R3',	'PIK3R4',	'PIK3R5'	,'PIK3R6',	'PIK3AP1',		'AKT1',	'AKT2',	'AKT3',		'IGF1',	'IGF1R',		'IL4',	'IL4R',		'NFATC1',	'NFATC2',	'NFATC3',	'NFATC4'	,'NFAT5',		'STAT6','CSF1','CSF1R')
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGA))
Marker_genes=intersect(Marker_genes,rownames(Gene_TCGA))

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

time1=Survival_CCGA[,4]
status1=Survival_CCGA[,5]
status1[is.na(status1)]=0
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.numeric(as.matrix(Gene_CCGA[Genes,])),dim(Gene_CCGA[Genes,]))

A=A4
# A4
# [1] -0.08057491  0.15705491  0.10197133 -0.03924833  6.38801401  0.27136837  0.14195919  0.01471263  0.39661697 -0.05663284  1.52427987
# [12]  0.02936325 -5.81479105 -0.06172801 -7.74822977 11.42249530 -0.23941925 -5.91291532 -1.91197970  0.04755584 -0.02362737

#### ROC for test set

# status1[status1==1&time1>365*5]=0  # 5-year survival

predicted1=matrix(0,1,dim(Gene_marker1)[2])
for (i in 1:dim(Gene_marker1)[2])
{
  predicted1[i]=sum(A*Gene_marker1[,i])
}

# predicted0=as.numeric(predicted1)
# predicted0=predicted1[!is.na(status1)]
# status0=na.omit(status1)
# pred <- prediction(as.numeric(predicted0),status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new = TRUE)
# plot(perf,colorize=FALSE, col="yellow") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )

# roc7=roc(status1,as.numeric(predicted1))
# AUC=auc(roc7)
# AUC


troc8<-timeROC(T=time1,delta=status1,
               marker=as.numeric(predicted1),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_4=predicted1


######### Age, Grade, Gender

## Read clinical data
setwd("Path")  
Survival_CGGA=read.csv("patient_survival data_CGGA.csv")

## Age
Age_CGGA=Survival_CGGA[,c(1,5)]
rownames(Age_CGGA)=Age_CGGA[,1]
Age=Age_CGGA[,-1]

Age[Age<=60]=1
Age[Age>60]=2
Age=as.matrix(t(Age))

# Gender
Gender_CGGA=Survival_CGGA[,c(1,4)]
rownames(Gender_CGGA)=Gender_CGGA[,1]
Gender=Gender_CGGA[,-1]
Gender=as.vector(Gender)

Gender[Gender=='F']=1
Gender[Gender=='M']=2
Gender=as.matrix(t(Gender))


# Grade

Grade_CGGA=Survival_CGGA[,c(1,7)]
rownames(Grade_CGGA)=Grade_CGGA[,1]
Grade=Grade_CGGA[,-1]
Grade=as.vector(Grade)

Grade[Grade==1|Grade==2]=1
Grade[Grade==3|Grade==4]=2
Grade=as.matrix(t(Grade))

################# Multi-variable COX model
Marker=rbind(as.numeric(Age),as.numeric(Gender),as.numeric(Grade),as.numeric(Score_1),as.numeric(t(Score_2)),as.numeric(Score_3),as.numeric(Score_4))

time1=Survival_CGGA[,8]
status1=Survival_CGGA[,9]
status1[is.na(status1)]=0
# status1[status1==1 & time1>365*5]=0  # 5-year survival

fit <- coxph(Surv(time1, status1) ~ t(Marker), method="breslow")

summary(fit)

################# Uni-variable COX model
time=time1
status=status1

Marker=rbind(as.numeric(Age))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Gender))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Grade))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)


Marker=rbind(Score_1)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(Score_2)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 

Marker=rbind(Score_3)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)  

Marker=rbind(Score_4)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)  

######## 5-year Survival evaluation (CGGA dataset)

status[status==1&time>365*5]=0  # 5-year survival status

################# Multi-variable COX model
Marker=rbind(as.numeric(Age),as.numeric(Gender),as.numeric(Grade),as.numeric(Score_1),as.numeric(t(Score_2)),as.numeric(Score_3),as.numeric(Score_4))

fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")

summary(fit)

################# Uni-variable COX model
Marker=rbind(as.numeric(Age))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Gender))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(as.numeric(Grade))
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)


Marker=rbind(Score_1)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)

Marker=rbind(Score_2)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 

Marker=rbind(Score_3)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)  

Marker=rbind(Score_4)
fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit) 

