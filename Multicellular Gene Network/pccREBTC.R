###PCC rebTC####
p_cor_rebTC<-array(data=NA,c(ncol(t_sel_rTC),ncol(t_sel_rTC),4))#??Ԥ???վ???
e_cor_rebTC<-array(data=NA,c(ncol(t_sel_rTC),ncol(t_sel_rTC),4))

for(k in 1:4)
{
  net_rebTC<-data.frame()
  for(i in 1:ncol(t_sel_rTC))#start cor.test single reb
  {for(j in 1:ncol(t_sel_rTC))
  {if(i<j)
  {cor<-cor.test(ver_sel_rTC[,i,k],ver_sel_rTC[,j,k])
  p_cor_rebTC[i,j,k]<-cor$p.value
  e_cor_rebTC[i,j,k]<-cor$estimate
  if((abs(e_cor_rebTC[i,j,k])>0.95)==TRUE&&(p_cor_rebTC[i,j,k]<0.05)==TRUE)
  {net_rebTC<-rbind(net_rebTC,data.frame(colnames(ver_sel_eTC)[i],colnames(ver_sel_eTC)[j],e_cor_rebTC[i,j,k],p_cor_rebTC[i,j,k]))}
  }}
    assign(paste("net_rebTC",k,sep=""),net_rebTC)
    write.csv(net_rebTC,file=paste("net_rebTC",k,".csv",sep=""))
    assign(paste("p_cor_rebTC",k,sep=""),p_cor_rebTC[,,k])
    assign(paste("e_cor_rebTC",k,sep=""),e_cor_rebTC[,,k])
    write.csv(p_cor_rebTC[,,k],file=paste("p_cor_rebTC",k,".csv",sep=""))
    write.csv(e_cor_rebTC[,,k],file=paste("e_cor_rebTC",k,".csv",sep=""))
  }}
