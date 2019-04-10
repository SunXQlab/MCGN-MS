###PCC rebTAM####
p_cor_rebTAM<-array(data=NA,c(ncol(t_sel_rTAM),ncol(t_sel_rTAM),4))#??Ԥ???վ???
e_cor_rebTAM<-array(data=NA,c(ncol(t_sel_rTAM),ncol(t_sel_rTAM),4))

for(k in 1:4)
{
  net_rebTAM<-data.frame()
  for(i in 1:ncol(t_sel_rTAM))#start cor.test single reb
  {for(j in 1:ncol(t_sel_rTAM))
  {if(i<j)
  {cor<-cor.test(ver_sel_rTAM[,i,k],ver_sel_rTAM[,j,k])
  p_cor_rebTAM[i,j,k]<-cor$p.value
  e_cor_rebTAM[i,j,k]<-cor$estimate
  if((abs(e_cor_rebTAM[i,j,k])>0.95)==TRUE&&(p_cor_rebTAM[i,j,k]<0.05)==TRUE)
  {net_rebTAM<-rbind(net_rebTAM,data.frame(colnames(ver_sel_eTAM)[i],colnames(ver_sel_eTAM)[j],e_cor_rebTAM[i,j,k],p_cor_rebTAM[i,j,k]))}
  }}
    assign(paste("net_rebTAM",k,sep=""),net_rebTAM)
    write.csv(net_rebTAM,file=paste("net_rebTAM",k,".csv",sep=""))
    assign(paste("p_cor_rebTAM",k,sep=""),p_cor_rebTAM[,,k])
    assign(paste("e_cor_rebTAM",k,sep=""),e_cor_rebTAM[,,k])
    write.csv(p_cor_rebTAM[,,k],file=paste("p_cor_rebTAM",k,".csv",sep=""))
    write.csv(e_cor_rebTAM[,,k],file=paste("e_cor_rebTAM",k,".csv",sep=""))
  }}
