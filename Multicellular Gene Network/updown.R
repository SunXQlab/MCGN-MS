up_down<-function(x,y)#x<-e_cor_rebTAM1,y<-e_cor_epTAM
{ud<-x
ud[which(abs(ud)-abs(y)>0)]<-1
ud[which(abs(ud)-abs(y)<0)]<-(-1)
return(ud)}

up_reb1<-up_down(e_cor_reb1,e_cor_ep)#TAMTC
up_reb2<-up_down(e_cor_reb2,e_cor_ep)
up_reb3<-up_down(e_cor_reb3,e_cor_ep)
up_reb4<-up_down(e_cor_reb4,e_cor_ep)
up_reb<-up_reb1+up_reb2+up_reb3+up_reb4

up_rebTAM1<-up_down(e_cor_rebTAM1,e_cor_epTAM)#interTAM
up_rebTAM2<-up_down(e_cor_rebTAM2,e_cor_epTAM)
up_rebTAM3<-up_down(e_cor_rebTAM3,e_cor_epTAM)
up_rebTAM4<-up_down(e_cor_rebTAM4,e_cor_epTAM)
up_rebTAM<-up_rebTAM1+up_rebTAM2+up_rebTAM3+up_rebTAM4

up_rebTC1<-up_down(e_cor_rebTC1,e_cor_epTC)#interTC
up_rebTC2<-up_down(e_cor_rebTC2,e_cor_epTC)
up_rebTC3<-up_down(e_cor_rebTC3,e_cor_epTC)
up_rebTC4<-up_down(e_cor_rebTC4,e_cor_epTC)
up_rebTC<-up_rebTC1+up_rebTC2+up_rebTC3+up_rebTC4

