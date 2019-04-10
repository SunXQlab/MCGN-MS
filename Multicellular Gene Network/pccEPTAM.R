#use ver_sel_eTAM eTAM.
net_epTAM<-data.frame()
p_cor_epTAM<-matrix(data=NA,nrow = ncol(ver_sel_eTAM),ncol = ncol(ver_sel_eTAM))#??Ԥ???վ???
e_cor_epTAM<-matrix(data=NA,nrow = ncol(ver_sel_eTAM),ncol = ncol(ver_sel_eTAM))

e<-0.95
p<-0.05
for(i in 1:ncol(ver_sel_eTAM))#start cor.test
{for(j in 1:ncol(ver_sel_eTAM))
{if(i<j)
{cor<-cor.test(ver_sel_eTAM[,i],ver_sel_eTAM[,j])
p_cor_epTAM[i,j]<-cor$p.value
e_cor_epTAM[i,j]<-cor$estimate
if((abs(e_cor_epTAM[i,j])>e)==TRUE&&(p_cor_epTAM[i,j]<p)==TRUE)
{net_epTAM<-rbind(net_epTAM,data.frame(colnames(ver_sel_eTAM)[i],colnames(ver_sel_eTAM)[j],e_cor_epTAM[i,j],p_cor_epTAM[i,j]))}
}}}
nrow(net_epTAM)

write.csv(p_cor_epTAM,file="p_cor_epTAM0.95_0.05.csv")
write.csv(e_cor_epTAM,file="e_cor_epTAM0.9_0.01.csv")
write.csv(net_epTAM,file="net_epTAM0.95_0.05.csv")