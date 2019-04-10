########### treatment response prediction: Targeted therapy  chemotherapy  Radiotherapy

###TCGA data_targeted therapy

setwd("Path")  

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


Survival_TCGA=Survival_TCGA[Survival_TCGA[,4]=='YES',]
dim(Survival_TCGA)

######################## Select genes for signature   #####################################
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG) 

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]


##################################################################################################################
###  3-year evaluation
##################################################################################################################

### Macrophage-related gene signature
A0=A
# A0=c(0.010721273, 0.003121212 -0.795150775, 0.001349858 -0.041226003, 0.016433866, 0.452024766, 1.144764257, 0.565811887, 0.072629383)
# A0=c(1.0107790, 1.0031261, 0.4515132, 1.0013508, 0.9596122, 1.0165696, 1.5714909, 3.1417006, 1.7608768, 1.0753319)
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]


########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival status

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A0*Gene_marker[,i])   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]

time_3=as.numeric(time)
status_3=as.numeric(status)
# time[time>365*3]=365*3
status_3[status_3==1&time_3>365*3]=0


status0=as.numeric(na.omit(status_3))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="red") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc0=roc(status0,predicted0)
AUC=roc0$auc
AUC

# troc0<-timeROC(T=time,delta=status,
#                marker=as.numeric(predicted),cause=1,
#                iid=TRUE,
#                times=c(3,5)*365,ROC=TRUE)

##########  EGFR
Marker_genes=c('EGFR')  

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

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
fit
cox.zph(fit)
A1=fit$coefficients

########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A1*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time_3=as.numeric(time)
status_3=as.numeric(status)
status_3[status_3==1&time_3>365*3]=0
status0=as.numeric(na.omit(status_3))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc1=roc(status0,predicted0)
AUC=auc(roc1)
AUC

# troc1<-timeROC(T=time,delta=status,
#                marker=as.numeric(predicted),cause=1,
#                iid=TRUE,
#                times=c(3,5)*365,ROC=TRUE)




##########  Cheng et al Neurology 2016; 86 (24)  
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112¨C118


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

# ### COX Model
# fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
# fit
# cox.zph(fit)
# A3=fit$coefficients

A3=c(-0.6718,0.1658, 0.2584, -0.1811, 0.1165, -0.4046, 0.1543, 0.1223)  # values in that paper

########  ROC for Training set
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A3*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time_3=as.numeric(time)
status_3=as.numeric(status)
status_3[status_3==1&time_3>365*3]=0
status0=as.numeric(na.omit(status_3))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="black") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc2=roc(status0,predicted0)
AUC=auc(roc2)
AUC

# troc2<-timeROC(T=time,delta=status,
#                marker=as.numeric(predicted),cause=1,
#                iid=TRUE,
#                times=c(3,5)*365,ROC=TRUE)

######### Compare two ROC
# install.packages("pROC")
# library("pROC")

# ci(roc0)
# ci(roc1)
# ci(roc2)
# ci(roc3)
# 
# roc.test(roc1,roc0,na.rm=TRUE,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs EGFR  p-value = 0.008582  AUC of roc1 AUC of roc2: 0.6475246   0.7802768
# roc.test(roc2,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value = 0.04875   AUC of roc1 AUC of roc2: 0.7155310   0.7802768


#### Quail et al. IGF1-PI3K pathways (Science 2016)

# Marker_genes=c('PI3K','AKT','IGF1','IGF1R','IL4','IL4R','NFAT5','STAT6')
Marker_genes=c('PIK3R1',	'PIK3R2'	,'PIK3R3',	'PIK3R4',	'PIK3R5'	,'PIK3R6',	'PIK3AP1',		'AKT1',	'AKT2',	'AKT3',		'IGF1',	'IGF1R',		'IL4',	'IL4R',		'NFATC1',	'NFATC2',	'NFATC3',	'NFATC4'	,'NFAT5',		'STAT6','CSF1','CSF1R')
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGA))
Marker_genes=intersect(Marker_genes,rownames(Gene_TCGA))



# ### Training dataset
# Survival_Gene=Gene_CGGA[t(Marker_genes),]
# colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
# 
# 
# Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
# Survival=as.matrix(Survival)
# Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
# status=as.numeric(Survival[,5])
# time=as.numeric(Survival[,4])
# Genes=rownames(na.omit(Survival_Gene))
# Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene)),dim(Survival_Gene))
# 
# 
# time=time[!is.na(status)]
# Gene_marker=Gene_marker[,!is.na(status)]
# status=status[!is.na(status)]

# ### COX Model
# fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
# fit
# cox.zph(fit)
# A4=fit$coefficients

A4=c(-0.08057491,0.15705491, 0.10197133, -0.03924833,  6.38801401,  0.27136837,  0.14195919,  0.01471263,  0.39661697, -0.05663284,  1.52427987,  0.02936325, -5.81479105, -0.06172801, -7.74822977, 11.42249530, -0.23941925, -5.91291532, -1.91197970,  0.04755584, -0.02362737)

# #### ROC for training set
# 
# predicted=matrix(0,1,dim(Gene_marker)[2])
# for (i in 1:dim(Gene_marker)[2])
# {
#   predicted[i]=sum(A4*Gene_marker[,i])
# }
# 
# predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# time=as.numeric(time)
# status=as.numeric(status)
# pred <- prediction(predicted0,status)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="green") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# roc3=roc(status,predicted0)
# AUC=roc3$auc
# AUC
# 

########  ROC for Testing set
Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

###
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene)),dim(Survival_Gene))
status[status==1&time>365*3]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A4*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time=as.numeric(time)
status=as.numeric(status[!is.na(status)])
pred <- prediction(predicted0,status)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="green") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc3=roc(status,predicted0)
AUC=roc3$auc
AUC

roc.test(roc1,roc0,na.rm=TRUE,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs EGFR  p-value = 0.008582  AUC of roc1 AUC of roc2: 0.6475246   0.7802768
roc.test(roc2,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value = 0.04875   AUC of roc1 AUC of roc2: 0.7155310   0.7802768
roc.test(roc3,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value = 0.04875   AUC of roc1 AUC of roc2: 0.7155310   0.7802768


##################################################################################################################
###  5-year evaluation
##################################################################################################################
setwd("Path")  

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

Survival_TCGA=Survival_TCGA[Survival_TCGA[,4]=='YES',]
dim(Survival_TCGA)

######################## Select genes for signature   #####################################
# MNB=c('PXDN' ,  'CDH6' ,  'SCN3A' , 'ANPEP' , 'SEMA6B', 'CCDC37' ,'PRRG1' , 'FANCA' , 'TMEM26' ,'NETO2' )

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]


### Macrophage-related gene signature

# A0=c(0.0111123294  0.0008619239 -0.8772969022  0.0016958263 -0.0423078649  0.0196739563  0.0013511639  0.8284922208  0.0026937362  1.1843625406   0.5722500646  0.1010323338)
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene))



time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]


########  ROC
# status[status==1&time>365*5]=0  # 5-year survival status

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A0*Gene_marker[,i])   # A0 is from macrophage-related prognostic signature, trained from TCGA data of all glioma patients
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]

time_5=as.numeric(time)
status_5=as.numeric(status)
status_5[status_5==1&time_5>365*5]=0


status0=as.numeric(na.omit(status_5))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="red") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc0=roc(status0,predicted0)
AUC=auc(roc0)
AUC

troc0<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

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

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
fit
cox.zph(fit)
A1=fit$coefficients

########  ROC 
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A1*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]

time_5=as.numeric(time)
status_5=as.numeric(status)
status_5[status_5==1&time_5>365*5]=0
status0=as.numeric(na.omit(status_5))

pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc1=roc(status0,predicted0)
AUC=auc(roc1)
AUC

troc1<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)




##########  Cheng et al Neurology 2016; 86 (24)  
Marker_genes=c('FOXO3', 'IL6', 'IL10', 'ZBTB16', 'CCL18', 'AIMP1', 'FCGR2B', 'MMP9')  # Bao et al., CNS Neuroscience & Therapeutics 20 (2014) 112¨C118


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

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
fit
cox.zph(fit)
A3=fit$coefficients

A3=c(-0.6718,0.1658, 0.2584, -0.1811, 0.1165, -0.4046, 0.1543, 0.1223)  # values in that paper

########  ROC 
# status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A3*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time_5=as.numeric(time)
status_5=as.numeric(status)
status_5[status_5==1&time_5>365*5]=0
status0=as.numeric(na.omit(status_5))
pred <- prediction(predicted0,status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="black") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc2=roc(status0,predicted0)
AUC=auc(roc2)
AUC

troc2<-timeROC(T=time,delta=status,
               marker=as.numeric(predicted),cause=1,
               iid=TRUE,
               times=c(3,5)*365,ROC=TRUE)

#  
# roc.test(roc1,roc0,na.rm=TRUE,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs EGFR   p-value = 0.03085  AUC of roc1 AUC of roc2: 0.6461574   0.7404375 
# roc.test(roc2,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value =0.03903   AUC of roc1 AUC of roc2: 0.6718330   0.7404375


#### Quail et al. IGF1-PI3K pathways (Science 2016)

Survival_Gene=Gene_TCGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

###
status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene)),dim(Survival_Gene))
status[status==1&time>365*5]=0  # 5-year survival

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A4*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time=as.numeric(time)
status=as.numeric(status[!is.na(status)])
pred <- prediction(predicted0,status)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
par(new = TRUE)
plot(perf,colorize=FALSE, col="green") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc3=roc(status,predicted0)
AUC=roc3$auc
AUC

roc.test(roc1,roc0,na.rm=TRUE,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs EGFR  p-value = 0.008582  AUC of roc1 AUC of roc2: 0.6475246   0.7802768
roc.test(roc2,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value = 0.04875   AUC of roc1 AUC of roc2: 0.7155310   0.7802768
roc.test(roc3,roc0,paired=TRUE,method="bootstrap",alternative="less",boot.n=1e3)  # MNB vs Cheng et al.  p-value = 0.04875   AUC of roc1 AUC of roc2: 0.7155310   0.7802768



# 
# dev.new()
# plot( troc0$FP[,1], troc0$TP[,1],type = "l", lty = 1, pch = NA, col = "red",ylab="", xlab="");
# lines( troc1$FP[,1], troc1$TP[,1],type = "l", lty = 1, pch = NA, col = "blue"); 
# # lines( troc2$FP[,1], troc2$TP[,1],type = "l", lty = 1, pch = NA, col = "green"); 
# lines( troc3$FP[,1], troc3$TP[,1],type = "l", lty = 1, pch = NA, col = "black"); 
# lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
# text.legend=c(paste("AUC of signature 1:",as.character(round(troc0$AUC[1],3))),paste("AUC of signature 2:",as.character(round(troc1$AUC[1],3))),paste("AUC of signature 3:",as.character(round(troc3$AUC[1],3))))
# legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("blue","red","black"),bty="n",ncol=1)
# title(main="3-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
# legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("red","blue","black"),bty="n",ncol=1)
# 
# dev.new()
# plot( troc0$FP[,2], troc0$TP[,2],type = "l", lty = 1, pch = NA, col = "red",ylab="", xlab="");
# lines( troc1$FP[,2], troc1$TP[,2],type = "l", lty = 1, pch = NA, col = "blue"); 
# # lines( troc2$FP[,2], troc2$TP[,2],type = "l", lty = 1, pch = NA, col = "green"); 
# lines( troc3$FP[,2], troc3$TP[,2],type = "l", lty = 1, pch = NA, col = "black"); 
# lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
# text.legend=c(paste("AUC of signature 1:",as.character(round(troc0$AUC[2],3))),paste("AUC of signature 2:",as.character(round(troc1$AUC[2],3))),paste("AUC of signature 3:",as.character(round(troc3$AUC[2],3))))
# legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("blue","red","black"),bty="n",ncol=1)
# title(main="5-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
# legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("red","blue","black"),bty="n",ncol=1)
# 
