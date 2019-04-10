
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


#### Quail et al. pathways (Science 2016)

# Marker_genes=c('PI3K','AKT','IGF1','IGF1R','IL4','IL4R','NFAT5','STAT6')
Marker_genes=c('PIK3R1',	'PIK3R2'	,'PIK3R3',	'PIK3R4',	'PIK3R5'	,'PIK3R6',	'PIK3AP1',		'AKT1',	'AKT2',	'AKT3',		'IGF1',	'IGF1R',		'IL4',	'IL4R',		'NFATC1',	'NFATC2',	'NFATC3',	'NFATC4'	,'NFAT5',		'STAT6','CSF1','CSF1R')
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGA))
Marker_genes=intersect(Marker_genes,rownames(Gene_TCGA))

Survival_Gene=Gene_CGGA[t(Marker_genes),]
colnames(Survival_Gene)=gsub(".", "-", colnames(Survival_Gene), fixed = TRUE)


Survival=Survival_CGGA[intersect(rownames(Survival_CGGA),colnames(Survival_Gene)),]
Survival=as.matrix(Survival)
Survival_Gene=Survival_Gene[,intersect(rownames(Survival),colnames(Survival_Gene))]

### Training dataset
status=as.numeric(Survival[,5])
time=as.numeric(Survival[,4])
Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.numeric(as.matrix(Survival_Gene)),dim(Survival_Gene))


time=time[!is.na(status)]
Gene_marker=Gene_marker[,!is.na(status)]
status=status[!is.na(status)]

### COX Model
fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
fit
cox.zph(fit)
A4=fit$coefficients

# A4
# [1] -0.08057491  0.15705491  0.10197133 -0.03924833  6.38801401  0.27136837  0.14195919  0.01471263  0.39661697 -0.05663284  1.52427987
# [12]  0.02936325 -5.81479105 -0.06172801 -7.74822977 11.42249530 -0.23941925 -5.91291532 -1.91197970  0.04755584 -0.02362737

#### ROC for trining set

predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A4*Gene_marker[,i])
}

predicted0=as.numeric(predicted)
predicted0=predicted0[!is.na(status)]
time=as.numeric(time)
status=as.numeric(status)
pred <- prediction(predicted0,status)
perf <- performance(pred,"tpr","fpr")  # calculate probabilities for TPR/FPR for predictions
performance(pred,"auc") # shows calculated AUC for model
dev.new()
plot(perf,colorize=FALSE, col="green") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

roc3=roc(status,predicted0)
AUC=roc3$auc
AUC


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
# status[status==1&time>365*5]=0  # 5-year survival

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
dev.new()
plot(perf,colorize=FALSE, col="green") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


roc3=roc(status,predicted0)
AUC=roc3$auc
AUC

