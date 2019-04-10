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
# DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
# DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
# DNG=t(as.matrix(DN_human[,1]))
# MNB=as.character(DNG)  

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### 
time=as.numeric(Survival[,3])
status=as.numeric(Survival[,2])


Genes=rownames(na.omit(Survival_Gene))
# Genes='ORC1'
Gene_marker= matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene[Genes,]))  #

time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

##### ROC for test set

# A=load('F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/Model_Coefficient.RData')
predicted1=matrix(0,1,dim(as.matrix(Gene_marker))[2])
for (i in 1:dim(as.matrix(Gene_marker))[2])
{
  predicted1[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc1=roc(as.character(status),as.numeric(predicted1))  #0.7595
# AUC=auc(roc1)
# AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CCGA datasets
groups=matrix(0,1,length(status))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(as.numeric(time)/365, as.numeric(status)) ~ groups, data = as.data.frame(Gene_marker))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))

ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = "nrisk_cumevents", risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Predicted",
           legend.labs = c("Sensitive","Resistant"))


##################################################### Time dependent ROC ##############################################
ROC.DSST<-timeROC(T=time,delta=status,
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
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red"),bty="n",ncol=1)


### 3-year survival prediction

time_3=as.numeric(time)
status_3=as.numeric(status)
# time[time>365*3]=365*3
status_3[status_3==1&time_3>365*3]=0

##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker))[2])
for (i in 1:dim(as.matrix(Gene_marker))[2])
{
  predicted1[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status_3)]
status0=na.omit(status_3)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc11=roc(as.character(status_3),as.numeric(predicted1))  #0.8676
# AUC=auc(roc1)
# AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc11$sensitivities,roc11$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CCGA datasets
groups=matrix(0,1,length(status))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(as.numeric(time_3)/365, as.numeric(status)) ~ groups, data = as.data.frame(Gene_marker))  # status for K_M analysis
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))


# ggsurvplot(fit,  size = 1,  # change line size
#            linetype =  c(1, 2), # change line type by groups
#            break.time.by = 50, # break time axis by 250
#            palette = c("#E7B800", "#2E9FDF"), # custom color palette
#            conf.int = TRUE, # Add confidence interval
#            pval = TRUE # Add p-value
# )

# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, data = as.data.frame(Gene_marker), linetype = "strata", pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = "nrisk_cumevents", risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Predicted",
           legend.labs = c("Sensitive","Resistant"))


### 5-year survival prediction
# time_5=as.numeric(Survival[,3])
# status_5=as.numeric(Survival[,2])
# Genes=rownames(na.omit(Survival_Gene))
# Gene_marker= matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene[Genes,]))  #
# 
# time_5=time_5[!is.na(status_5)]
# Gene_marker=Gene_marker[,!is.na(status_5)]
# status_5=status_5[!is.na(status_5)]
# 
# time_5=as.numeric(time_5)
# status_5=as.numeric(status_5)
# # time[time>365*3]=365*3

time_5=as.numeric(time)
status_5=as.numeric(status)
status_5[status_5==1&time_5>365*5]=0

##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker))[2])
for (i in 1:dim(as.matrix(Gene_marker))[2])
{
  predicted1[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status_5)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
# performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc11=roc(as.character(status_5),as.numeric(predicted1))  #0.8676
# AUC=auc(roc1)
# AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CCGA datasets
groups=matrix(0,1,length(status))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(as.numeric(time_5)/365, as.numeric(status)) ~ groups, data = as.data.frame(Gene_marker))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))


# ggsurvplot(fit,  size = 1,  # change line size
#            linetype =  c(1, 2), # change line type by groups
#            break.time.by = 50, # break time axis by 250
#            palette = c("#E7B800", "#2E9FDF"), # custom color palette
#            conf.int = TRUE, # Add confidence interval
#            pval = TRUE # Add p-value
# )

# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, data = as.data.frame(Gene_marker), linetype = "strata", pval = TRUE, conf.int = TRUE,
           xlab='Time (Years)',
           risk.table = "nrisk_cumevents", risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Predicted",
           legend.labs = c("Sensitive","Resistant"))









################ TCGA chemotherapy ###########################

Gene_GBM=read.csv("TCGA_GBM_GeneExpression.csv")
rownames(Gene_GBM)=Gene_GBM[,1]
colnames(Gene_GBM)=gsub(".", "-", colnames(Gene_GBM), fixed = TRUE)
Gene_GBM=Gene_GBM[,-1]

Survival_GBM=read.csv("TCGA_GBM_Survival data.csv")
Survival_GBM=Survival_GBM[,c(1,33,32,62)]
rownames(Survival_GBM)=Survival_GBM[,1]
colnames(Survival_GBM)=c('SampleID','OS_event','OS_time','ChemoTherapy')

Gene_TCGA=Gene_GBM
colnames(Gene_TCGA)=gsub(".", "-", colnames(Gene_TCGA), fixed = TRUE)
Survival_TCGA=Survival_GBM

Survival_TCGA=Survival_TCGA[Survival_TCGA[,4]=='YES',]
dim(Survival_TCGA)


######################## Select genes for signature   #####################################
DN=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Table 1.csv")  # Differential gene network
DN_human=read.csv("F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Human Genes in Differential Network.csv")  # Differential gene network
DNG=t(as.matrix(DN_human[,1]))
MNB=as.character(DNG)  # Macrophage-associated gene signature
# [1] "WFDC1"   "PTK7"    "PTGER2"  "PXDN"    "MIP"     "TANC1"   "CDH6"    "SCN3A"   "LTC4S"   "ANPEP"   "SEMA6B"  "ABCG4"   "GPR171"  "REP15"   "XCR1"    "SPSB4"  
# [17] "COL19A1" "EPS8L2"  "BTLA"    "CCDC37"  "PTPRF"   "DPP4"    "ARG1"    "PRRG1"   "GPNMB"   "FANCA"   "TMEM26"  "SH3TC2"  "NETO2"

# ID=ensemble2symbol[intersect(rownames(ensemble2symbol),GeneName),] 

Survival_Gene=Gene_TCGA[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_TCGA[intersect(rownames(Survival_TCGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### 
time=as.numeric(Survival[,3])
status=Survival[,2]


Genes=rownames(na.omit(Survival_Gene))

Gene_marker= matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene[Genes,]))  #


##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker))[2])
for (i in 1:dim(as.matrix(Gene_marker))[2])
{
  predicted1[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc11=roc(as.character(status),as.numeric(predicted1))  #0.7595
AUC=auc(roc1)
AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CCGA datasets
groups=matrix(0,1,length(status))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(as.numeric(time)/365, as.numeric(status)) ~ groups, data = as.data.frame(Gene_marker))
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


##################################################### Time dependent ROC ##############################################
ROC.DSST<-timeROC(T=time,delta=as.numeric(status),
                  marker=as.numeric(predicted1),cause=1,
                  weighting="cox",
                  times=c(1,3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[1],3))),paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[2],3))),paste("AUC at 7 years:",as.character(round(ROC.DSST$AUC[3],3))))
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)

################# LDH1 Wild type/mutation subtype of Glioma patients

setwd("Path")  

Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]

Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,67)]
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','LDH1')

Survival_LGG=as.matrix(Survival_LGG)
Survival_LDH1wd=Survival_LGG[Survival_LGG[,4]=='NO',]
Survival_LDH1mt=Survival_LGG[Survival_LGG[,4]=='YES',]

Gene_LGG=as.matrix(Gene_LGG)
Gene_LDH1wd=Gene_LGG[Survival_LGG[,4]=='NO',]
Gene_LDH1mt=Gene_LGG[Survival_LGG[,4]=='YES',]

Survival_Gene=Gene_LDH1wd[t(MNB),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_LDH1wd[intersect(rownames(Survival_LDH1wd),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### 
time=as.numeric(Survival[,3])
status=Survival[,2]


Genes=rownames(na.omit(Survival_Gene))
# Genes='ORC1'
Gene_marker= matrix(as.numeric(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene[Genes,]))  #


##### ROC for test set
predicted1=matrix(0,1,dim(as.matrix(Gene_marker))[2])
for (i in 1:dim(as.matrix(Gene_marker))[2])
{
  predicted1[i]=sum(A*Gene_marker[,i])
}

predicted0=as.numeric(predicted1)
predicted0=predicted1[!is.na(status)]
status0=na.omit(status)
pred <- prediction(as.numeric(predicted0),status0)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="blue") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc11=roc(as.character(status),as.numeric(predicted1))  #0.7595
# AUC=auc(roc1)
# AUC
## optimal combination
opt <- which.max(rowSums(cbind(roc1$sensitivities,roc1$specificities)))
## optimal cut-off point 
sort(predicted1,F)[opt]

##########  K-M survival curves for MNB in TCGA and CCGA datasets
groups=matrix(0,1,length(status))
groups[predicted1<=sort(predicted1,F)[opt]]=1
groups[predicted1>sort(predicted1,F)[opt]]=2
# groups[predicted0<=median(predicted0)]=1
# groups[predicted0>median(predicted0)]=2
groups=t(groups)
groups=as.numeric(groups)


fit<- survfit(Surv(as.numeric(time)/365, as.numeric(status)) ~ groups, data = as.data.frame(Gene_marker))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))

ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.y.text.col = TRUE)

