####PCC####
#use sel_t_reb_TAM etc.

ver_sel_eTAM<-t(sel100_eTAM)
ver_sel_eTC<-t(sel100_eTC)

####PCC epTAMTC####
#use ver_sel_eTAM etc.
net_ep<-data.frame()
p_cor_ep<-matrix(data=NA,nrow = ncol(ver_sel_eTAM),ncol = ncol(ver_sel_eTC))#??Ԥ???վ???
e_cor_ep<-matrix(data=NA,nrow = ncol(ver_sel_eTAM),ncol = ncol(ver_sel_eTC))

e<-0.95
p<-0.05
for(i in 1:ncol(ver_sel_eTAM))#start cor.test
{for(j in 1:ncol(ver_sel_eTC))
{cor<-cor.test(ver_sel_eTAM[,i],ver_sel_eTC[,j])
p_cor_ep[i,j]<-cor$p.value
e_cor_ep[i,j]<-cor$estimate
if((abs(e_cor_ep[i,j])>e)==TRUE&&(p_cor_ep[i,j]<p)==TRUE)
{net_ep<-rbind(net_ep,data.frame(colnames(ver_sel_eTAM)[i],colnames(ver_sel_eTC)[j],e_cor_ep[i,j],p_cor_ep[i,j]))}
}}
nrow(net_ep)

write.csv(p_cor_ep,file="p_cor_ep0.9_0.05.csv")
write.csv(e_cor_ep,file="e_cor_ep0.9_0.05.csv")
write.csv(net_ep,file="net_ep0.95_0.05.csv")


