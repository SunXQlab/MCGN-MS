#use ver_sel_eTC eTC.
net_epTC<-data.frame()
p_cor_epTC<-matrix(data=NA,nrow = ncol(ver_sel_eTC),ncol = ncol(ver_sel_eTC))#??Ԥ???վ???
e_cor_epTC<-matrix(data=NA,nrow = ncol(ver_sel_eTC),ncol = ncol(ver_sel_eTC))

e<-0.95
p<-0.05
for(i in 1:ncol(ver_sel_eTC))#start cor.test
{for(j in 1:ncol(ver_sel_eTC))
{if(i<j)
{cor<-cor.test(ver_sel_eTC[,i],ver_sel_eTC[,j])
p_cor_epTC[i,j]<-cor$p.value
e_cor_epTC[i,j]<-cor$estimate
if((abs(e_cor_epTC[i,j])>e)==TRUE&&(p_cor_epTC[i,j]<p)==TRUE)
{net_epTC<-rbind(net_epTC,data.frame(colnames(ver_sel_eTC)[i],colnames(ver_sel_eTC)[j],e_cor_epTC[i,j],p_cor_epTC[i,j]))}
}}}
nrow(net_epTC)

write.csv(p_cor_epTC,file="p_cor_epTC0.9_0.01.csv")
write.csv(e_cor_epTC,file="e_cor_epTC0.9_0.01.csv")
write.csv(net_epTC,file="net_epTC0.95_0.05.csv")