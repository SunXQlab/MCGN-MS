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
Gene_CGGA=Gene_CGGA[,!is.na(status)]
status=status[!is.na(status)]



# ########  ROC for Training set
# predicted=matrix(0,1,dim(Gene_marker)[2])
# for (i in 1:dim(Gene_marker)[2])
# {
#   predicted[i]=sum(A*Gene_marker[,i])
# }
# 
# predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# status0=as.matrix(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# dev.new()
# plot(perf,colorize=FALSE, col="red") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# roc0=roc((status0),predicted0)
# roc0
# AUC=auc(roc0)
# AUC


####### Random 12 gene set 
dev.new()
G=matrix(as.numeric(Gene_CGGA),dim(Gene_CGGA))
Gene_set=Gene_CGGA[apply(G,1,function(x) all(x)>=1e-3),]
AUC=matrix(0,1,1000)
for (j in 1:1000)
{
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
  Gene_CGGA=Gene_CGGA[,!is.na(status)]
  status=status[!is.na(status)]
  
  
Gene_marker_random=Gene_CGGA[sample(1:dim(Gene_CGGA)[1],12),]
Gene_name=rownames(Gene_marker_random)
Gene_marker_random=matrix(as.double(as.matrix(Gene_marker_random)),dim(Gene_marker_random))
fit <- coxph(Surv(time, status) ~ t(Gene_marker_random), method="breslow")
A_random=coefficients(fit)
########  ROC for Training set
# predicted=matrix(0,1,dim(Gene_marker_random)[2])
# for (i in 1:dim(Gene_marker_random)[2])
# {
#   predicted[i]=sum(as.numeric(A_random)*Gene_marker_random[,i])
# }
# 
# predicted0=as.numeric(predicted)
# predicted0=predicted0[!is.na(status)]
# status0=as.matrix(na.omit(status))
# pred <- prediction(predicted0,status0)
# perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
# par(new=TRUE)
# plot(perf,colorize=FALSE, col="blue") # plot ROC curve
# lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# 
# roc0=roc((status0),as.numeric(predicted0))
# roc0
# AUC[j]=auc(roc0)

### Test datasets
Survival_Gene=Gene_TCGA[t(Gene_name),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker1=matrix(as.double(as.matrix(Survival_Gene[Gene_name,])),dim(Survival_Gene[Gene_name,]))



time1=time[!is.na(status)]
status1=status[!is.na(status)]
Genes=rownames(na.omit(Survival_Gene))

Gene_marker1=Gene_marker1[,!is.na(status)]

Gene_marker1[is.na(Gene_marker1)]=rexp(sum(is.na(Gene_marker1)), rate=1)*100

predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
for (i in 1:dim(as.matrix(Gene_marker1))[2])
{
  predicted1[i]=sum(A_random*Gene_marker1[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status1)]
status0=na.omit(status1)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
if (j%%10==0)
{
par(new=TRUE)
plot(perf,colorize=FALSE, col="#7F7FFF", yaxt="n",xaxt="n",xlab="",ylab="",axes=FALSE) # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
}
roc1=roc(status1,as.numeric(predicted0))
AUC[j]=auc(roc1)

print(j)
}
AUC


########## 

####### Random 12 gene set from DEG with FC >1.5
DEG_name=read.csv('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEGs_fc_over_1.5.csv')
# MNB=read.csv('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEGs.csv')

DEG_name=(as.matrix(DEG_name[,3]))
DEG_name=DEG_name[DEG_name!=""]
DEG=intersect(rownames(Gene_CGGA),t(DEG_name))

DEG=intersect(DEG,rownames(Gene_CGGA))
DEG=intersect(DEG,rownames(Gene_TCGA))

G=matrix(as.numeric(Gene_CGGA),dim(Gene_CGGA))
Gene_set=Gene_CGGA[apply(G,1,function(x) all(x)>=1e-3),]
AUC=matrix(0,1,1000)

# dev.new()
for (j in 1:1000)
{
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
  Gene_DEG=Gene_CGGA[DEG,!is.na(status)]
  status=status[!is.na(status)]
  
  
  Gene_marker_random=Gene_DEG[sample(1:dim(Gene_DEG)[1],12),]
  Gene_name=rownames(Gene_marker_random)
  Gene_marker_random=matrix(as.double(as.matrix(Gene_marker_random)),dim(Gene_marker_random))
  fit <- coxph(Surv(time, status) ~ t(Gene_marker_random), method="breslow")
  A_random=coefficients(fit)
  ########  ROC for Training set
  # predicted=matrix(0,1,dim(Gene_marker_random)[2])
  # for (i in 1:dim(Gene_marker_random)[2])
  # {
  #   predicted[i]=sum(as.numeric(A_random)*Gene_marker_random[,i])
  # }
  # 
  # predicted0=as.numeric(predicted)
  # predicted0=predicted0[!is.na(status)]
  # status0=as.matrix(na.omit(status))
  # pred <- prediction(predicted0,status0)
  # perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  # performance(pred,"auc") # shows calculated AUC for model
  # par(new=TRUE)
  # plot(perf,colorize=FALSE, col="blue") # plot ROC curve
  # lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  # 
  # roc0=roc((status0),as.numeric(predicted0))
  # roc0
  # AUC[j]=auc(roc0)
  
  ### Test datasets
  Survival_Gene=Gene_TCGA[t(Gene_name),]
  colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)
  
  
  Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
  Survival=as.matrix(Survival)
  Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]
  
  status=as.numeric(Survival[,2])
  time=as.numeric(Survival[,3])
  Genes=rownames(na.omit(Survival_Gene))
  Gene_marker1=matrix(as.double(as.matrix(Survival_Gene[Gene_name,])),dim(Survival_Gene[Gene_name,]))
  
  
  
  time1=time[!is.na(status)]
  status1=status[!is.na(status)]
  Genes=rownames(na.omit(Survival_Gene))
  
  Gene_marker1=Gene_marker1[,!is.na(status)]
  
  Gene_marker1[is.na(Gene_marker1)]=rexp(sum(is.na(Gene_marker1)), rate=1)*100
  
  predicted1=matrix(0,1,dim(as.matrix(Gene_marker1))[2])
  for (i in 1:dim(as.matrix(Gene_marker1))[2])
  {
    predicted1[i]=sum(A_random*Gene_marker1[,i])
  }
  
  predicted0=as.numeric(predicted1)
  predicted0=predicted1[!is.na(status1)]
  status0=na.omit(status1)
  pred <- prediction(as.numeric(predicted0),status0)
  perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
  performance(pred,"auc") # shows calculated AUC for model
  if (j%%10==0)
  {
    par(new=TRUE)
    plot(perf,colorize=FALSE, col="#7F7FFF", yaxt="n",xaxt="n",xlab="",ylab="",axes=FALSE) # plot ROC curve
    lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  }
  roc1=roc(status1,as.numeric(predicted0))
  AUC[j]=roc1$auc
  
  print(j)
}
AUC


##################################################################################################################
##################################################################################################################

### macrophage-related gene signature 
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG)  # Macrophage-associated gene signature

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
par(new=TRUE)
plot(perf,colorize=FALSE, col="#FF0000", lwd=3 ) # plot ROC curve
lines(c(0,1),c(0,1),col = "black", lty = 4 )

roc1=roc(status1,as.numeric(predicted0))
AUC0=roc1$auc
AUC0

dev.new()
hist(AUC)
par(new=TRUE)
plot(AUC0,0:10)

dev.new()
hist(AUC, 20, freq = FALSE)
curve(dnorm(x, mean = mean(AUC), sd = sd(AUC)), min(AUC), max(AUC), add = TRUE, col = "blue", lwd = 2 , main='')


# Kernel Density Plot


d <- density(AUC) # returns the density data
dev.new()
plot(d, main="Prob. distribution",xlab="AUC values") # plots the result

x=AUC
h<-hist(x, breaks=15, col="#7F7FFF", xlab="AUC values",
        main="")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

dev.new()
hist(AUC,breaks=100,prob=TRUE, col="#7F7FFF", xlab="AUC values", main="",xlim = c(0.4, 0.9),ylim = c(0, 10),xaxs = "i", yaxs ="i")
lines(density(AUC), col="darkblue", lwd=2)
par(new=TRUE)
lines(rep(AUC0,11),0:10,col="#FF0000", lwd=3)

p.value<-length(AUC[AUC>AUC0])/1000  # 1e-3

p.value
                          
                          
