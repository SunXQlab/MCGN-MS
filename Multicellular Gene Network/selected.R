# setwd("F:/Path")


# DATA<-read.csv("DATA02.csv")#read file  
# rownames(DATA)<-DATA[,1]
# DATA<-DATA[,-1]
# head(DATA)
# 
# DATA<-as.matrix(DATA)
# write.csv(DATA,file="Gene_expression_DATA.csv")
# save(DATA, file="F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Submission/Journal of Translational Medicine/Code/Gene_expression_DATA.RData")

# DATA has been saved in R workshop and can be loaded to directly use. 

seqrebTAM<-seq(from=1,to=7,by=2)#read datalist
reb_TAM<-DATA[,seqrebTAM]
seqrebTC<-seq(from=2,to=8,by=2)
reb_TC<-DATA[,seqrebTC]
seqepTAM<-seq(from=9,to=20,by=2)
ep_TAM<-DATA[,seqepTAM]
seqepTC<-seq(from=10,to=20,by=2)
ep_TC<-DATA[,seqepTC]
seqvehTAM<-seq(from=21,to=30,by=2)
veh_TAM<-DATA[,seqvehTAM]
seqvehTC<-seq(from=22,to=31,by=2)
veh_TC<-DATA[,seqvehTC]

mean_rTAM<-apply(reb_TAM,1,mean)
mean_rTC<-apply(reb_TC,1,mean)
mean_eTAM<-apply(ep_TAM,1,mean)
mean_eTC<-apply(ep_TC,1,mean)
mean_vTAM<-apply(veh_TAM,1,mean)
mean_vTC<-apply(veh_TC,1,mean)

log_rTAM<-log2(reb_TAM)/log2(mean_vTAM)
log_eTAM<-log2(ep_TAM)/log2(mean_vTAM)
log_rTC<-log2(reb_TC)/log2(mean_vTC)
log_eTC<-log2(ep_TC)/log2(mean_vTC)

#####t test TAM########################
t_reb_TAM<-log_rTAM#create t var
t_ep_TAM<-log_eTAM

P_value=matrix(NA,1,nrow(t_reb_TAM))
i_TAM<-c() 
for (i in 1:nrow(t_reb_TAM))#start t_TAM
{if(all(t_reb_TAM[i,]==t_reb_TAM[i,1])==FALSE&&
    all(t_ep_TAM[i,]==t_ep_TAM[i,1])==FALSE&&
    is.finite(t_reb_TAM[i,1])&&is.finite(t_reb_TAM[i,2])&&is.finite(t_reb_TAM[i,3])&&is.finite(t_reb_TAM[i,4])&&
    is.finite(t_ep_TAM[i,1])&&is.finite(t_ep_TAM[i,2])&&is.finite(t_ep_TAM[i,3])&&is.finite(t_ep_TAM[i,4])&&is.finite(t_ep_TAM[i,5])&&is.finite(t_ep_TAM[i,6]))
{P_value[i]<-t.test(t_reb_TAM[i,],t_ep_TAM[i,])$p.value
if(P_value[i]<0.05)
{i_TAM<-c(i_TAM,i)
}
}
}

P_adjust=p.adjust(P_value[i_TAM],method="fdr",n=length(P_value[i_TAM]))
i_TAM=i_TAM[P_adjust<0.05]
length(i_TAM)


# A=which(!is.na(P_value))
# P_adjust=p.adjust(P_value[A],method="fdr",n=length(A))
# i_TAM=A[P_adjust<0.05]
# length(i_TAM)


sel_t_rTAM<-t_reb_TAM[i_TAM,]
sel_t_eTAM<-t_ep_TAM[i_TAM,]

######fold change select100TAM#########
mean_rTAM<-apply(sel_t_rTAM,1,mean)
mean_eTAM<-apply(sel_t_eTAM,1,mean)

fc_TAM<-abs(log2(mean_rTAM/mean_eTAM))
order_fc_TAM<-order(fc_TAM,decreasing = TRUE)[1:50]
sel100_rTAM<-sel_t_rTAM[order_fc_TAM,]
sel100_eTAM<-sel_t_eTAM[order_fc_TAM,]
rownames(sel100_eTAM)

# #write.csv(i_TAM,file="t_i_TAM100.csv")#save i
# write.csv(sel100_rTAM,file="selected50rebTAM.csv")
# write.csv(sel100_eTAM,file="selected050epTAM.csv")
#######################



#####t test TAM########################
t_reb_TC<-log_rTC#create t var
t_ep_TC<-log_eTC

P_value=matrix(0,1,nrow(t_reb_TAM))
i_TC<-c() 
#pv<-c()
for (i in 1:nrow(t_reb_TC))#start t_TC
{if(all(t_reb_TC[i,]==t_reb_TC[i,1])==FALSE&&
    all(t_ep_TC[i,]==t_ep_TC[i,1])==FALSE&&
    is.finite(t_reb_TC[i,1])&&is.finite(t_reb_TC[i,2])&&is.finite(t_reb_TC[i,3])&&is.finite(t_reb_TC[i,4])&&
    is.finite(t_ep_TC[i,1])&&is.finite(t_ep_TC[i,2])&&is.finite(t_ep_TC[i,3])&&is.finite(t_ep_TC[i,4])&&is.finite(t_ep_TC[i,5])&&is.finite(t_ep_TC[i,6]))
{P_value[i]<-t.test(t_reb_TC[i,],t_ep_TC[i,])$p.value
if(P_value[i]<0.05)
{i_TC<-c(i_TC,i)
}}
}

P_adjust=p.adjust(P_value[i_TC],method="fdr",n=length(P_value[i_TC]))
i_TC=i_TC[P_adjust<0.05]
length(i_TC)

sel_t_rTC<-t_reb_TC[i_TC,]
sel_t_eTC<-t_ep_TC[i_TC,]

######fold change select100TC#########
mean_rTC<-apply(sel_t_rTC,1,mean)
mean_eTC<-apply(sel_t_eTC,1,mean)

fc_TC<-abs(log2(mean_rTC/mean_eTC))


DEG=unique(c(names(fc_TC[fc_TC>log2(1.5)]),names(fc_TAM[fc_TAM>log(1.5)])))  
write.csv(DEG,'F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEGs_fc_over_1.5.csv')

# DEG=cbind(ep_TC[intersect(DEG,rownames(ep_TC)),],ep_TAM[intersect(DEG,rownames(ep_TAM)),],reb_TAM[intersect(DEG,rownames(reb_TAM)),],reb_TC[intersect(DEG,rownames(reb_TC)),])
DEG_data=DATA[intersect(rownames(DATA),rownames(DEG)),]
write.csv(DEG_data,'F:/Glioma Drug Resistance Gene Network/Glioma Drug Resistance Gene Network/Survival analysis - New/DEG_expression_30mouse.csv')


order_fc_TC<-order(fc_TC,decreasing = TRUE)[1:50]
sel100_rTC<-sel_t_rTC[order_fc_TC,]
sel100_eTC<-sel_t_eTC[order_fc_TC,]

# write.csv(sel100_rTC,file="selected50rebTC.csv")
# write.csv(sel100_eTC,file="selected50epTC.csv")
#######################
###heatmap###
library(ggplot2)
heatmap(cbind(sel100_eTAM,sel100_rTAM))
heatmap(cbind(sel100_eTC,sel100_rTC))
