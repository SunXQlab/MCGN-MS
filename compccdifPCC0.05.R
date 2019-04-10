compcc<-function(xe,ye,xp,yp)#xe<-e_cor_reb1,ye<-e_cor_ep,xp<-p_cor_reb1,yp<-p_cor_ep
{a<-abs(xe-ye)>0.05
#which(x==TRUE,arr.ind = TRUE)
b<-xp<0.05
c<-yp<0.05
d<-a*b*c
return(d)
}


compcc_reb1<-compcc(e_cor_reb1,e_cor_ep,p_cor_reb1,p_cor_ep)
compcc_reb2<-compcc(e_cor_reb2,e_cor_ep,p_cor_reb2,p_cor_ep)
compcc_reb3<-compcc(e_cor_reb3,e_cor_ep,p_cor_reb3,p_cor_ep)
compcc_reb4<-compcc(e_cor_reb4,e_cor_ep,p_cor_reb4,p_cor_ep)
compcc_reb<-compcc_reb1*compcc_reb2*compcc_reb3*compcc_reb4
n_compcc_reb<-which(compcc_reb==1,arr.ind = TRUE)
net_compcc_reb<-data.frame(paste("TAM",rownames(sel100_eTAM)[n_compcc_reb[,1]]),paste("TC",rownames(sel100_eTC)[n_compcc_reb[,2]]),up_reb[n_compcc_reb])
write.csv(net_compcc_reb,file="net_compcc_reb_pcc0.05.csv")


compcc_rebTAM1<-compcc(e_cor_rebTAM1,e_cor_epTAM,p_cor_rebTAM1,p_cor_epTAM)
compcc_rebTAM2<-compcc(e_cor_rebTAM2,e_cor_epTAM,p_cor_rebTAM2,p_cor_epTAM)
compcc_rebTAM3<-compcc(e_cor_rebTAM3,e_cor_epTAM,p_cor_rebTAM3,p_cor_epTAM)
compcc_rebTAM4<-compcc(e_cor_rebTAM4,e_cor_epTAM,p_cor_rebTAM4,p_cor_epTAM)
compcc_rebTAM<-compcc_rebTAM1*compcc_rebTAM2*compcc_rebTAM3*compcc_rebTAM4
n_compcc_rebTAM<-which(compcc_rebTAM==1,arr.ind = TRUE)
net_compcc_rebTAM<-data.frame(paste("TAM",rownames(sel100_eTAM)[n_compcc_rebTAM[,1]]),paste("TAM",rownames(sel100_eTAM)[n_compcc_rebTAM[,2]]),up_rebTAM[n_compcc_rebTAM])
write.csv(net_compcc_rebTAM,file="net_compcc_rebTAM_pcc0.05.csv")

compcc_rebTC1<-compcc(e_cor_rebTC1,e_cor_epTC,p_cor_rebTC1,p_cor_epTC)
compcc_rebTC2<-compcc(e_cor_rebTC2,e_cor_epTC,p_cor_rebTC2,p_cor_epTC)
compcc_rebTC3<-compcc(e_cor_rebTC3,e_cor_epTC,p_cor_rebTC3,p_cor_epTC)
compcc_rebTC4<-compcc(e_cor_rebTC4,e_cor_epTC,p_cor_rebTC4,p_cor_epTC)
compcc_rebTC<-compcc_rebTC1*compcc_rebTC2*compcc_rebTC3*compcc_rebTC4
n_compcc_rebTC<-which(compcc_rebTC==1,arr.ind = TRUE)
net_compcc_rebTC<-data.frame(paste("TC",rownames(sel100_eTC)[n_compcc_rebTC[,1]]),paste("TC",rownames(sel100_eTC)[n_compcc_rebTC[,2]]),up_rebTC[n_compcc_rebTC])
write.csv(net_compcc_rebTC,file="net_compcc_rebTC_pcc0.05.csv")
