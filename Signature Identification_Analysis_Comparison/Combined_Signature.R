########### treatment response prediction: Targeted therapy  chemotherapy  Radiotherapy


# DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
# DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
# DNG=t(as.matrix(DN_human[,1]))
# MNB=as.character(DNG)  # Macrophage-associated gene signature



##################################################################################################################
### Training datasets  (CGGA dataset)
##################################################################################################################

#### Age
setwd("F:/2017学生论文/scRNA-seq Intercellular signaling pathway")  
Survival_CGGA=read.csv("patient_survival data_CGGA.csv")

## Age
Age_CGGA=Survival_CGGA[,c(1,5)]
rownames(Age_CGGA)=Age_CGGA[,1]
Age=Age_CGGA[,-1]

predicted=matrix(0,1,length(Age))
for (i in 1:length(Age))
{
  predicted[i]=Age[i]   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

time=Survival_CGGA[,8]
status=Survival_CGGA[,9]
status[is.na(status)]=0

troc1<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

#### Grade
Grade_CGGA=Survival_CGGA[,c(1,7)]
rownames(Grade_CGGA)=Grade_CGGA[,1]
Grade=Grade_CGGA[,-1]
Grade=as.vector(Grade)

Grade[Grade==1|Grade==2]=1
Grade[Grade==3|Grade==4]=2
Grade=as.matrix(t(Grade))
# fit <- coxph(Surv(time, status) ~ (as.numeric(Grade)), method="breslow")
# A=coef(fit)
A= 2.251227
predicted=matrix(0,1,length(Grade))
for (i in 1:length(Grade))
{
  predicted[i]=sum(A*as.numeric(Grade[,i]))   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

troc2<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)


### Signature 1 (Macrophage-related genes)
Marker_genes=c("PXDN" ,  "CDH6", "SCN3A", "ANPEP", "SEMA6B", "CCDC37", "DPP4",  "PRRG1", "GPNMB", "FANCA", "TMEM26", "NETO2" )

Survival_Gene=Gene_CCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_CCGA[intersect(rownames(Survival_CCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,5])
time=as.numeric(Survival[,4])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)

# time=time[!is.na(status)]
# Gene_marker=Gene_marker[,!is.na(status)]
# status=status[!is.na(status)]


A0=c(0.0111123294, 0.0008619239, -0.8772969022, 0.0016958263, -0.0423078649, 0.0196739563, 0.0013511639, 0.8284922208, 0.0026937362, 1.1843625406, 0.5722500646, 0.1010323338)

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A0*Gene_marker[,i])   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}


troc3<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_1=predicted

# ##########  Signature 3 ### Cheng et al. 
################################################################
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112–118

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

troc4<-timeROC(T=time1,delta=status1,
               marker=as.numeric(predicted1),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_3=predicted1

##########  Combined Signature #############


################ LASSO regression for multivariate COX model
Marker=rbind(as.numeric(Age),as.numeric(Grade),Score_1,Score_3)

# surv=Surv(as.double(time),as.double(status))
# cvfit<-cv.glmnet(t(Marker),surv,family="cox",alpha=1)
# plot(cvfit)
# A=coef(cvfit, s = "lambda.min")
# A=coef(cvfit, s="lambda.1se")
# which(A!=0)
# sum(A!=0)
# # A_combined=A[A!=0]
# A=A_combined

fit <- coxph(Surv(time, status) ~ t(Marker), method="breslow")
summary(fit)
A=coefficients(fit)

########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Marker)[2])
for (i in 1:dim(Marker)[2])
{
  predicted[i]=sum(A*Marker[,i])
}


troc5<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

######### Compare two ROC

compare(troc5,troc1)  # ompute p-values of comparison test
compare(troc5,troc2)  # ompute p-values of comparison test
compare(troc5,troc3)  # ompute p-values of comparison test
compare(troc5,troc4)  # ompute p-values of comparison test

dev.new()
plot( troc1$FP[,1], troc1$TP[,1],type = "l", lty = 1, pch = NA, col = "brown",ylab="", xlab="");
lines( troc2$FP[,1], troc2$TP[,1],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,1], troc3$TP[,1],type = "l", lty = 1, pch = NA, col = "MediumBlue"); 
lines( troc4$FP[,1], troc4$TP[,1],type = "l", lty = 1, pch = NA, col = "black"); 
lines( troc5$FP[,1], troc5$TP[,1],type = "l", lty = 1, pch = NA, col = "red");
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of Age:",as.character(round(troc1$AUC[1],3))),paste("AUC of Grade:",as.character(round(troc2$AUC[1],3))),paste("AUC of signature 1:",as.character(round(troc3$AUC[1],3))),paste("AUC of signature 3:",as.character(round(troc4$AUC[1],3))),paste("AUC of combined signature:",as.character(round(troc5$AUC[1],3))))
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)
title(main="3-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)

dev.new()
plot( troc1$FP[,2], troc1$TP[,2],type = "l", lty = 1, pch = NA, col = "brown",ylab="", xlab="");
lines( troc2$FP[,2], troc2$TP[,2],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,2], troc3$TP[,2],type = "l", lty = 1, pch = NA, col = "MediumBlue"); 
lines( troc4$FP[,2], troc4$TP[,2],type = "l", lty = 1, pch = NA, col = "black"); 
lines( troc5$FP[,2], troc5$TP[,2],type = "l", lty = 1, pch = NA, col = "red");
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of Age:",as.character(round(troc1$AUC[2],3))),paste("AUC of Grade:",as.character(round(troc2$AUC[2],3))),paste("AUC of signature 1:",as.character(round(troc3$AUC[2],3))),paste("AUC of signature 3:",as.character(round(troc4$AUC[2],3))),paste("AUC of combined signature:",as.character(round(troc5$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)
title(main="5-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)


########################  Test TCGA dataset

setwd("F:/2017学生论文/scRNA-seq Intercellular signaling pathway")  

Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,93)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Target_Therapy')

#####
Gene_GBM=read.csv("TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,109)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','Target_Therapy')

Gene_TCGA=cbind(Gene_GBM[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),],Gene_LGG[intersect(rownames(Gene_GBM),rownames(Gene_LGG)),])
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=rbind(Survival_GBM,Survival_LGG)

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]


##################################################################################################################
##################################################################################################################

### Training dataset
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))
colnames(Gene_marker)=colnames(Survival_Gene)



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]


##################################################################################################################
##################################################################################################################
#### Age
## Age
setwd("F:/2017学生论文/scRNA-seq Intercellular signaling pathway")  
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")

Age_LGG=Survival_LGG[,c(1,20)]
rownames(Age_LGG)=Age_LGG[,1]

Age_GBM=Survival_GBM[,c(1,56)]
rownames(Age_GBM)=Age_GBM[,1]                     


Age_TCGA=c(Age_GBM[,2],Age_LGG[,2])
names(Age_TCGA)=c(rownames(Age_GBM),rownames(Age_LGG))
Age=Age_TCGA[PatientNo_1]

predicted=matrix(0,1,length(Age))
for (i in 1:length(Age))
{
  predicted[i]=Age[i]   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

troc1<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

#### Grade
Grade_TCGA=c(rep('GBM',length(rownames(Gender_GBM))),rep('LGG',length(rownames(Gender_LGG))))
names(Grade_TCGA)=c(rownames(Gender_GBM),rownames(Gender_LGG))
Grade=Grade_TCGA[PatientNo_1]

Grade[Grade=='LGG']=1
Grade[Grade=='GBM']=2

Grade=as.matrix(t(Grade))

fit <- coxph(Surv(time, status) ~ (as.numeric(Grade)), method="breslow")
A=coef(fit)
predicted=matrix(0,1,length(Grade))
for (i in 1:length(Grade))
{
  predicted[i]=sum(A*as.numeric(Grade[,i]))   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

troc2<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)


### Signature 1 (Macrophage-related genes)
Marker_genes=c("PXDN" ,  "CDH6", "SCN3A", "ANPEP", "SEMA6B", "CCDC37", "DPP4",  "PRRG1", "GPNMB", "FANCA", "TMEM26", "NETO2" )

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


A0=c(0.0111123294, 0.0008619239, -0.8772969022, 0.0016958263, -0.0423078649, 0.0196739563, 0.0013511639, 0.8284922208, 0.0026937362, 1.1843625406, 0.5722500646, 0.1010323338)
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A0*Gene_marker[,i])   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}



troc3<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)
Score_1=predicted

##########  Cheng et al Neurology 2016; 86 (24)  
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112–118


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


troc4<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

Score_3=predicted

##########  Combined Signature #############

################ LASSO regression for multivariate COX model
Marker=rbind(as.numeric(Age), as.numeric(Grade), Score_1, Score_3)

predicted=matrix(0,1,dim(Marker)[2])
for (i in 1:dim(Marker)[2])
{
  predicted[i]=sum(A*Marker[,i])
}


troc5<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

######### Compare two ROC

compare(troc5,troc1)  # ompute p-values of comparison test
compare(troc5,troc2)  # ompute p-values of comparison test
compare(troc5,troc3)  # ompute p-values of comparison test
compare(troc5,troc4)  # ompute p-values of comparison test


dev.new()
plot( troc1$FP[,1], troc1$TP[,1],type = "l", lty = 1, pch = NA, col = "brown",ylab="", xlab="");
lines( troc2$FP[,1], troc2$TP[,1],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,1], troc3$TP[,1],type = "l", lty = 1, pch = NA, col = "MediumBlue"); 
lines( troc4$FP[,1], troc4$TP[,1],type = "l", lty = 1, pch = NA, col = "black"); 
lines( troc5$FP[,1], troc5$TP[,1],type = "l", lty = 1, pch = NA, col = "red");
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of Age:",as.character(round(troc1$AUC[1],3))),paste("AUC of Grade:",as.character(round(troc2$AUC[1],3))),paste("AUC of signature 1:",as.character(round(troc3$AUC[1],3))),paste("AUC of signature 3:",as.character(round(troc4$AUC[1],3))),paste("AUC of combined signature:",as.character(round(troc5$AUC[1],3))))
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)
title(main="3-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)

dev.new()
plot( troc1$FP[,2], troc1$TP[,2],type = "l", lty = 1, pch = NA, col = "brown",ylab="", xlab="");
lines( troc2$FP[,2], troc2$TP[,2],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,2], troc3$TP[,2],type = "l", lty = 1, pch = NA, col = "MediumBlue"); 
lines( troc4$FP[,2], troc4$TP[,2],type = "l", lty = 1, pch = NA, col = "black"); 
lines( troc5$FP[,2], troc5$TP[,2],type = "l", lty = 1, pch = NA, col = "red");
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of Age:",as.character(round(troc1$AUC[2],3))),paste("AUC of Grade:",as.character(round(troc2$AUC[2],3))),paste("AUC of signature 1:",as.character(round(troc3$AUC[2],3))),paste("AUC of signature 3:",as.character(round(troc4$AUC[2],3))),paste("AUC of combined signature:",as.character(round(troc5$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)
title(main="5-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,1,1),legend=text.legend,col=c("brown","green","MediumBlue","black","red"),bty="n",ncol=1)




