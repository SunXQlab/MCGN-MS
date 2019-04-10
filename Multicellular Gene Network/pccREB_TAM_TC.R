####PCC single reb/ep####
#ver_sel_eTAM<-t(sel_eTAM)#vert
#ver_sel_eTC<-t(sel_eTC)

####set array####
t_sel_rTAM<-t(sel100_rTAM)
t_sel_rTC<-t(sel100_rTC)
#t_sel_eTAM<-t(sel_t_ep_TAM)
#t_sel_eTC<-t(sel_t_ep_TC)

ver_sel_rTAM<-array(NA,c(7,ncol(t_sel_rTAM),4))#create matrix adding single reb
ver_sel_rTC<-array(NA,c(7,ncol(t_sel_rTC),4))
for(i in 1:4)
{ver_sel_rTAM[,,i]<-rbind(ver_sel_eTAM,t_sel_rTAM[i,])
ver_sel_rTC[,,i]<-rbind(ver_sel_eTC,t_sel_rTC[i,])
}

###beginning:rebPCC####
p_cor_reb<-array(data=NA,c(ncol(t_sel_rTAM),ncol(t_sel_rTC),4))#??Ԥ???վ???
e_cor_reb<-array(data=NA,c(ncol(t_sel_rTAM),ncol(t_sel_rTC),4))

for(k in 1:4)
{
  net_reb<-data.frame()
  for(i in 1:ncol(t_sel_rTAM))#start cor.test single reb
  {for(j in 1:ncol(t_sel_rTC))
  {cor<-cor.test(ver_sel_rTAM[,i,k],ver_sel_rTC[,j,k])
  p_cor_reb[i,j,k]<-cor$p.value
  e_cor_reb[i,j,k]<-cor$estimate
  if((abs(e_cor_reb[i,j,k])>0.95)==TRUE&&(p_cor_reb[i,j,k]<0.05)==TRUE)
  {net_reb<-rbind(net_reb,data.frame(colnames(ver_sel_eTAM)[i],colnames(ver_sel_eTC)[j],e_cor_reb[i,j,k],p_cor_reb[i,j,k]))}
  }}
  assign(paste("net_reb",k,sep=""),net_reb)
  write.csv(net_reb,file=paste("net_reb",k,".csv",sep=""))
  assign(paste("p_cor_reb",k,sep=""),p_cor_reb[,,k])
  assign(paste("e_cor_reb",k,sep=""),e_cor_reb[,,k])
  write.csv(p_cor_reb[,,k],file=paste("p_cor_reb",k,".csv",sep=""))
  write.csv(e_cor_reb[,,k],file=paste("e_cor_reb",k,".csv",sep=""))
}
