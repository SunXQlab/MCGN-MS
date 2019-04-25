
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

# install.packages("rms")
library(rms)


detach("package:org.Hs.eg.db")
detach("package:foreign")

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
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG)  # Macrophage-associated gene signature
# [1] "WFDC1"   "PTK7"    "PTGER2"  "PXDN"    "MIP"     "TANC1"   "CDH6"    "SCN3A"   "LTC4S"   "ANPEP"   "SEMA6B"  "ABCG4"   "GPR171"  "REP15"   "XCR1"    "SPSB4"  
# [17] "COL19A1" "EPS8L2"  "BTLA"    "CCDC37"  "PTPRF"   "DPP4"    "ARG1"    "PRRG1"   "GPNMB"   "FANCA"   "TMEM26"  "SH3TC2"  "NETO2"

# ID=ensemble2symbol[intersect(rownames(ensemble2symbol),GeneName),] 


##################################################################################################################
##################################################################################################################

### Training dataset
Survival_Gene=Gene_CGGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,5])
time=as.numeric(Survival[,4])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

##########  rbsurv ##############
# source("http://bioconductor.org/biocLite.R")
# biocLite("rbsurv")
# library(rbsurv)

surv=Surv(as.double(time),as.double(status))

cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1,nfolds=10)
# dev.new()
# plot(cvfit)
A=coef(cvfit, s = "lambda.min")
# A=coef(cvfit, s="lambda.1se")
A=as.numeric(A)
A[(A!=0)]
which(A!=0)
sum(A!=0)
MNB[which(A!=0)]

# save(cvfit,file='F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/cvfit.RData')
# load('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/cvfit.RData')

# 
# save(A,file='F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/Model_Coefficient.RData')
# load('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/Model_Coefficient.RData')


#########  LASSO 



# DEG=unique(c(names(fc_TC[fc_TC>log2(1.5)]),names(fc_TAM[fc_TAM>log(1.5)])))  
# write.csv(DEG,'F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEGs_fc_over_1.5.csv')

DEG=read.csv('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEGs_fc_over_1.5.csv')

DEG=(as.matrix(DEG[,3]))
DEG=DEG[DEG!=""]
DEG=intersect(rownames(Gene_CGGA),t(DEG))

DEG=intersect(DEG,rownames(Gene_CGGA))
DEG=intersect(DEG,rownames(Gene_TCGA))

### Training dataset
Survival_Gene=Gene_CGGA[t(DEG),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,5])
time=as.numeric(Survival[,4])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

surv=Surv(as.double(time),as.double(status))

cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1,nfolds=10)
# dev.new()
# plot(cvfit)
A1=coef(cvfit, s = "lambda.min")
# A=coef(cvfit, s="lambda.1se")
A1=as.numeric(A1)
A1[(A1!=0)]
which(A1!=0)
sum(A1!=0)
DEG[which(A1!=0)]
A1

save(cvfit,file='F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/cvfit_DEG_LASSO.RData')
# load('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/cvfit.RData')
 
save(A1,file='F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/Model_Coefficient_DEG_LASSO.RData')
# load('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/Model_Coefficient.RData')


########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A1*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
status0=as.matrix(na.omit(status))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc0=roc((status0),predicted0)
roc0
AUC=roc0$auc
AUC
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
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))

########### time dependent ROC 
ROC.DSST<-timeROC(T=time,delta=status,
                  marker=as.numeric(predicted),cause=1,
                  weighting="cox",
                  times=c(3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)



##################################################################################################################
##################################################################################################################

### Test datasets
Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time1=time[!is.na(status)]
status1=status[!is.na(status)]
Genes=rownames(na.omit(Survival_Gene))

Gene_marker1=Gene_marker1[,!is.na(status)]
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

roc1=roc(status1,as.numeric(predicted0))
AUC=roc1$auc
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


# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))

####### Time dependent ROC
ROC.DSST<-timeROC(T=time1,delta=status1,
                  marker=as.numeric(predicted1),cause=1,
                  weighting="cox",
                  times=c(3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
# lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)

############### Robusteness 

## time
## status
## predicted

AUC_OS=matrix(0,1,100)
AUC_3Y=matrix(0,1,100)
AUC_5Y=matrix(0,1,100)

for (i in 1:100)
{
  
  SampleID=sample(1:length(status1),floor(length(status1)*0.5))
  time_2=time1[SampleID]
  status_2=status1[SampleID]
  predicted_2=predicted0[SampleID]
  
  roc_2=roc(status_2,predicted_2)
  AUC=roc_2$auc
  AUC_OS[i]=AUC
  
  
  troc_2<-timeROC(T=time_2,delta=status_2,
                  marker=as.numeric(predicted_2),cause=1,
                  iid=TRUE,
                  times=c(3,5)*365,ROC=TRUE)
  
  
  troc_2$AUC
  
  AUC_3Y[i]=troc_2$AUC[1]
  AUC_5Y[i]=troc_2$AUC[2]
  
}

dev.new()
hist(AUC_OS)
hist(AUC_3Y)
hist(AUC_5Y)

###### Test for LASSO

Survival_Gene=Gene_TCGA[t(DEG),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
# Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.double(as.matrix(Survival_Gene)),dim(Survival_Gene))



time1=time[!is.na(status)]
status1=status[!is.na(status)]
# Genes=rownames(na.omit(Survival_Gene))

Gene_marker1=Gene_marker1[,!is.na(status)]
Gene_marker1[is.na(Gene_marker1)]=0


##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A1*Gene_marker1[,i])
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

roc1=roc(status1,as.numeric(predicted0))
AUC=roc1$auc
AUC


ROC.DSST<-timeROC(T=time1,delta=status1,
                  marker=as.numeric(predicted0),cause=1,
                  weighting="cox",
                  times=c(3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
# lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)



#########################

##### Robustness
AUC_OS_LASSO=matrix(0,1,100)
AUC_3Y_LASSO=matrix(0,1,100)
AUC_5Y_LASSO=matrix(0,1,100)

for (i in 1:100)
{
  
  SampleID=sample(1:length(status1),floor(length(status1)*0.5))
  time_2=time1[SampleID]
  status_2=status1[SampleID]
  predicted_2=predicted0[SampleID]
  
  roc_rand=roc(status_2,predicted_2)
  AUC=as.numeric(roc_rand$auc)
  AUC_OS_LASSO[i]=AUC
  
  
  troc_2<-timeROC(T=time_2,delta=status_2,
                  marker=as.numeric(predicted_2),cause=1,
                  iid=TRUE,
                  times=c(3,5)*365,ROC=TRUE)
  
  
  troc_2$AUC
  
  AUC_3Y_LASSO[i]=troc_2$AUC[1]
  AUC_5Y_LASSO[i]=troc_2$AUC[2]
  
}

# hist(AUC_OS)
# hist(AUC_3Y)
# hist(AUC_5Y)

#### Wilcoxon Test computing p value

wilcox.test(AUC_OS, AUC_OS_LASSO, alternative = "greater", paired=T) 
wilcox.test(AUC_3Y, AUC_3Y_LASSO, alternative = "greater", paired=T) 
wilcox.test(AUC_5Y, AUC_5Y_LASSO, alternative = "greater", paired=T) 


# overall survival ROC AUC
dev.new()
hist(AUC_OS_LASSO,breaks=10,prob=TRUE, col="LightGreen", border = "gray", xlab="AUC values", main="",xlim = c(0.65, 0.9),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_OS,breaks=10,prob=TRUE, col = "#7F7FFF",  border = "gray", xlab="AUC values", main="",xlim = c(0.65, 0.9),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_OS_LASSO), col="#4D9221", lwd=3)
lines(density(AUC_OS), col="blue", lwd=3)


# 3-year survival ROC AUC
dev.new()
hist(AUC_3Y_LASSO,breaks=10,prob=TRUE, col="LightGreen", border = "gray", xlab="AUC values", main="",xlim = c(0.7, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_3Y,breaks=10,prob=TRUE, col="#7F7FFF", border = "gray", xlab="AUC values", main="",xlim = c(0.7, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_3Y_LASSO), col="#4D9221", lwd=3)
lines(density(AUC_3Y), col="blue", lwd=3)

# 5-year survival ROC AUC
dev.new()
hist(AUC_5Y_LASSO,breaks=10,prob=TRUE, col="LightGreen", border = "gray", xlab="AUC values", main="",xlim = c(0.7, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_5Y,breaks=10,prob=TRUE, col="#7F7FFF", border = "gray", xlab="AUC values", main="",xlim = c(0.7, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_5Y_LASSO), col="#4D9221", lwd=3)
lines(density(AUC_5Y), col="blue", lwd=3)



