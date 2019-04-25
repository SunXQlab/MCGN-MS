
###### Install packages in R 
# source("http://bioc.ism.ac.jp/biocLite.R")
# biocLite("AnnotationDbi")
# install.packages("AnnotationDbi")
# install.packages("bit")
# biocLite("org.Hs.eg.db")
# library(bit)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# ###### Survival package
# install.packages("lattice")
# library(lattice)
# 
# 
# ######## plot KM curves
# install.packages("reshape2")
# library(reshape2)
# install.packages("data.table")
# library(data.table)
# install.packages("zoo")
# library("zoo")
# install.packages("survminer")
# library("survminer")
# install.packages("survival")
# library("survival")
# 
# #### ROC package
# install.packages("gplots")
# library(gplots)
# install.packages("ROCR")
# library(ROCR)
# install.packages("dplyr")
# library("dplyr")
# install.packages("lubridate")
# library(lubridate)
# install.packages("robustbase")
# install.packages("caret")
# library(caret)
# install.packages("pROC")
# library(pROC)
# 
# install.packages("timeROC")
# install.packages("prodlim")
# library(prodlim)
# library(quantreg)
# install.packages("quantreg")
# install.packages("polspline")
# library("polspline")
# library(timeROC)
# install.packages("foreign")
# library(foreign)
# 
# # install.packages("rms")
# library(rms)
# 
# 
# detach("package:org.Hs.eg.db")
# detach("package:foreign")
# 
# install.packages("glmnet")
# library(glmnet)  # manually load this package by downloading and installing from Toll. 


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

######################## DEGs   #####################################
DEG_data=read.csv('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEG_expression_30mouse.csv')
rownames(DEG_data)=DEG_data[,1]
DEG_data=DEG_data[,-1]
DEG=rownames(DEG_data)

#### Map mouse gene symbol to human's
setwd('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Submission/Journal of Translational Medicine/Code')
H2M=read.delim("HOM_MouseHumanSequence.rpt")
mouse=H2M[H2M[,2]=="mouse, laboratory",c(1,4)]
human=H2M[H2M[,2]=="human",c(1,4)]
H2M=merge(mouse,human,by.x=1,by.y=1)[,2:3]
H2M=H2M[!duplicated(H2M[,1])&!duplicated(H2M[,2]),]
H2M=as.matrix(H2M)
rownames(H2M)=H2M[,1]
DEG=intersect(DEG,H2M[,1])

DEG_h2m=H2M[DEG,]

DEG_hum=DEG_data[rownames(DEG_h2m),]
rownames(DEG_hum)=DEG_h2m[,2]


dim(DEG_hum)  #1141   30
######################## Correlation network   #####################################
X=DEG_hum
PC= matrix(NA,nrow=dim(X)[1],dim(X)[1])
PC_p=matrix(NA,nrow=dim(X)[1],dim(X)[1])

for (i in 511:dim(X)[1])
{
  print(i)
  for (j in 1:dim(X)[1])
  {
  B=cor.test(as.numeric(X[i,]),as.numeric(X[j,]),method="pearson")
  PC[i,j]=B$estimate
  PC_p[i,j]=B$p.value    
  }
}
PC_saved=PC

PC[PC_p>=0.05]=0

S=matrix(0,1,dim(X)[1])
for (i in 1:dim(X)[1])
{
  S[i]=sum(abs(PC[i,]))
}

dev.new(); hist(S) 


Node_selected=rownames(DEG_hum[order(S,decreasing = TRUE)[1:100],])
Node_selected=intersect(Node_selected,rownames(Gene_CGGA))
Node_selected=intersect(Node_selected,rownames(Gene_TCGA))
Node_selected=Node_selected[1:12]

######################## Select genes for signature   #####################################
MNB=Node_selected

### Training dataset
Survival_Gene=Gene_CGGA[(MNB),]
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

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
A_PCCNet=coefficients(fit)

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A_PCCNet*Gene_marker[,i])
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
Gene_marker1=matrix(as.double(as.matrix(Survival_Gene)),dim(Survival_Gene))



time1=time[!is.na(status)]
status1=status[!is.na(status)]
Genes=rownames(na.omit(Survival_Gene))

Gene_marker1=Gene_marker1[,!is.na(status)]
##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A_PCCNet*Gene_marker1[,i])
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

#### Macrophage-related gene signature
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG)  # Macrophage-associated gene signature

# A
# [1]  0.0000000000  0.0000000000  0.0000000000  0.0111123294  0.0000000000  0.0000000000  0.0008619239 -0.8772969022  0.0000000000  0.0016958263
# [11] -0.0423078649  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0196739563
# [21]  0.0000000000  0.0013511639  0.0000000000  0.8284922208  0.0026937362  1.1843625406  0.5722500646  0.0000000000  0.1010323338

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

###
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

###### Test for PCCNet (Correlation Network based on DEGS)
PCCNet=Node_selected
Survival_Gene=Gene_TCGA[t(PCCNet),]
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
  predicted1[i]=sum(A_PCCNet*Gene_marker1[,i])
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
AUC_OS_PCCNet=matrix(0,1,100)
AUC_3Y_PCCNet=matrix(0,1,100)
AUC_5Y_PCCNet=matrix(0,1,100)

for (i in 1:100)
{
  
  SampleID=sample(1:length(status1),floor(length(status1)*0.5))
  time_2=time1[SampleID]
  status_2=status1[SampleID]
  predicted_2=predicted0[SampleID]
  
  roc_rand=roc(status_2,predicted_2)
  AUC=as.numeric(roc_rand$auc)
  AUC_OS_PCCNet[i]=AUC
  
  
  troc_2<-timeROC(T=time_2,delta=status_2,
                  marker=as.numeric(predicted_2),cause=1,
                  iid=TRUE,
                  times=c(3,5)*365,ROC=TRUE)
  
  
  troc_2$AUC
  
  AUC_3Y_PCCNet[i]=troc_2$AUC[1]
  AUC_5Y_PCCNet[i]=troc_2$AUC[2]
  
}

# hist(AUC_OS)
# hist(AUC_3Y)
# hist(AUC_5Y)

#### WilcoxonTest computing p value

wilcox.test(AUC_OS, AUC_OS_PCCNet, alternative = "greater") 
wilcox.test(AUC_3Y, AUC_3Y_PCCNet, alternative = "greater") 
wilcox.test(AUC_5Y, AUC_5Y_PCCNet, alternative = "greater") 


# overall survival ROC AUC
dev.new()
hist(AUC_OS_PCCNet,breaks=10,prob=TRUE, col="#B8E186", border = "gray", xlab="AUC values", main="",xlim = c(0.5, 0.9),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_OS,breaks=10,prob=TRUE, col = "#7F7FFF",  border = "gray", xlab="AUC values", main="",xlim = c(0.5, 0.9),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_OS_PCCNet), col="#4D9221", lwd=3)
lines(density(AUC_OS), col="blue", lwd=3)


# 3-year survival ROC AUC
dev.new()
hist(AUC_3Y_PCCNet,breaks=10,prob=TRUE, col="#B8E186", border = "gray", xlab="AUC values", main="",xlim = c(0.5, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_3Y,breaks=10,prob=TRUE, col="#7F7FFF", border = "gray", xlab="AUC values", main="",xlim = c(0.5, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_3Y_PCCNet), col="#4D9221", lwd=3)
lines(density(AUC_3Y), col="blue", lwd=3)

# 5-year survival ROC AUC
dev.new()
hist(AUC_5Y_PCCNet,breaks=10,prob=TRUE, col="#B8E186", border = "gray", xlab="AUC values", main="",xlim = c(0.5, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
par(new=TRUE)
hist(AUC_5Y,breaks=10,prob=TRUE, col="#7F7FFF", border = "gray", xlab="AUC values", main="",xlim = c(0.5, 1),ylim = c(0, 30),xaxs = "i", yaxs ="i")
lines(density(AUC_5Y_PCCNet), col="#4D9221", lwd=3)
lines(density(AUC_5Y), col="blue", lwd=3)



