# 


# Load --------------------------------------------------------------------
sdate=as.Date("2004-05-01");# Begining of Florida WY2005
edate=as.Date("2018-04-30");# End of Florida WY2018

#tmap_mode("view")
tm_shape(enp.slough.clip)+tm_polygons(alpha=0.5)+
  tm_shape(sfwmd.mon.q)+tm_dots(size=0.05)+
  tm_shape(fce.sites)+tm_dots(col="yellow")+
  tm_shape(sfwmd.mon.wq)+tm_dots(col="red")
#tm_shape(sfwmd.mon.wx)+tm_dots(col="blue",size=0.07)

HW.flow.sites=data.frame(SITE=c(paste0("S12",LETTERS[1:4]),"S333","S334","S332D","S199","S18C","S197"),
                         DBKEY=c("01313","00610","00621","01310","15042","FB752","TA413","88897","15760","15763"),
                         Region=c(rep("SRS",6),rep("TSPh",3),"Tide"))
HW.wq.sites=data.frame(SITE=c(paste0("S12",LETTERS[1:4]),"S333",rep("S332D",2),"S199","S18C","S197"),
                       WQSite=c(paste0("S12",LETTERS[1:4]),"S333","S332D","S332DX","S177","S18C","S197"))
flow.wq.xwalk=rbind(HW.wq.sites,data.frame(SITE=c("S333R","S18C_R"),WQSite=c("S333","S18C")))
wq.params=data.frame(param=c("TP","TN","NOx","TKN","DOC","TOC"),Test.Number=c(25,80,18,21,89,100))

HW.wq=data.frame()
for(i in 1:nrow(HW.wq.sites)){
  tmp=SFWMD.DBHYDRO.Data.WQ(sdate,edate,HW.wq.sites$WQSite[i],wq.params$Test.Number)
  HW.wq=rbind(HW.wq,tmp)
  print(i)
}
unique(HW.wq$Station.ID)
HW.wq=subset(HW.wq,Collection.Method=="G")
HW.wq$AbsMDL=with(HW.wq,ifelse(Value<0,abs(Value),Value))
HW.wq=merge(HW.wq,wq.params,"Test.Number")
HW.wq=merge(HW.wq,HW.wq.sites,by.x="Station.ID",by.y="WQSite")

HW.wq.xtab=cast(HW.wq,Station.ID+SITE+Date.EST~param,value="AbsMDL",mean)
HW.wq.xtab$TN.final=with(HW.wq.xtab,SFWMD.TN.Combine(NOx,TKN,TN))

HW.q=data.frame()
for(i in 1:nrow(HW.flow.sites)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,HW.flow.sites$DBKEY[i])
  tmp$DBKEY=as.character(HW.flow.sites$DBKEY[i])
  HW.q=rbind(HW.q,tmp)
  print(i)
}
HW.q=merge(HW.q,HW.flow.sites,by=c("DBKEY"))
HW.q$Data.Value=with(HW.q,ifelse(Data.Value<0,0,Data.Value))
HW.q$Date.EST=date.fun(HW.q$Date)
HW.q$WY=WY(HW.q$Date.EST)

HW.flow.da.xtab=cast(HW.q,Date.EST+WY~SITE,value="Data.Value",mean)
HW.flow.da.xtab$S333R=with(HW.flow.da.xtab,ifelse(S333-S334<0,0,S333-S334))
HW.flow.da.xtab$S18C_R=with(HW.flow.da.xtab,ifelse(S18C-S197<0,0,S18C-S197))
HW.flow.da.xtab$tflow.srs=rowSums(HW.flow.da.xtab[,c(paste0("S12",LETTERS[1:4]),"S333R")],na.rm=T)
HW.flow.da.xtab$tflow.TSv1=rowSums(HW.flow.da.xtab[,c("S332D","S199","S18C_R")],na.rm=T)
HW.flow.da.xtab$tflow.TSv2=rowSums(HW.flow.da.xtab[,c("S332D","S199","S18C")],na.rm=T)

plot(tflow.srs~Date.EST,subset(HW.flow.da.xtab,WY==2018),type="l")
plot(S18C~Date.EST,subset(HW.flow.da.xtab,WY==2018),type="l")
plot(S197~Date.EST,subset(HW.flow.da.xtab,WY==2018),type="l")

#xlim.val=date.fun(c("2004-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years");
xlim.val=date.fun(c("2017-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 months");xmin=seq(xlim.val[1],xlim.val[2],"1 months");
ylim.val=c(0,7e6);by.y=1e6;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"TSPh_flow.tiff"),width=6.5,height=4.75,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,2.5,0.25,0.25),oma=c(1,1.5,1,0.5),mgp=c(3,1,0));
layout(matrix(c(1:2),2,1,byrow=T),heights=c(1,0.2))

plot(tflow.TSv2~Date.EST,HW.flow.da.xtab,type="n",ylim=ylim.val,xlim=xlim.val,ylab=NA,xlab=NA,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(HW.flow.da.xtab,shaded.range(Date.EST,rep(0,length(Date.EST)),cfs.to.m3d(tflow.TSv2),"grey50",lty=1))
with(HW.flow.da.xtab,shaded.range(Date.EST,rep(0,length(Date.EST)),cfs.to.m3d(tflow.TSv1),"indianred1",lty=1))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),1);axis_fun(2,ymaj,ymin,ymaj/1e5,1);box(lwd=1)
mtext(side=2,line=2,"Daily Discharge (x 10\u2075; m\u00B3 Yr\u207B\u00B9)")
mtext(side=1,line=1.75,"Date (Month-Year)")
mtext(side=3,"TS/Ph")
abline(v=c(irma.landfall.FL),col="black",lty=2)

plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
legend.text=c("S332D+S199+S18C","S332D+S199+(S18C-S197)")
pt.cols=c(adjustcolor("grey50",0.25),adjustcolor("indianred1",0.25))
legend(0.5,0.5,legend=legend.text,lty=c(0,0),col=c("grey50","indianred1"),lwd=c(0.5,0.5),pch =c(22,22),pt.bg=pt.cols,pt.cex=1.25,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

#with(HW.flow.da.xtab,lines(Date.EST,tflow.TS,col="Red"))
#abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="black",lty=2)

vars=c("S12A","S12B","S12C","S12D","S333R","S18C_R","S199","S332D","S197")
HW.flow.da.melt=melt(data.frame(HW.flow.da.xtab[,c("Date.EST","WY",vars)]),id.vars=c("Date.EST","WY"),variable_name="SITE")
HW.flow.da.melt=merge(HW.flow.da.melt,flow.wq.xwalk,"SITE")
unique(HW.flow.da.melt$SITE)

HW.flow.wq.da=merge(HW.flow.da.melt,data.frame(HW.wq.xtab),by.x=c("Date.EST","SITE"),by.y=c("Date.EST","Station.ID"),all.x=T)
HW.flow.wq.da=HW.flow.wq.da[order(HW.flow.wq.da$SITE,HW.flow.wq.da$Date.EST),]
HW.flow.wq.da$TP.inter=with(HW.flow.wq.da,ave(TP,SITE,FUN=function(x)dat.interp(x)))
HW.flow.wq.da$TN.inter=with(HW.flow.wq.da,ave(TN.final,SITE,FUN=function(x)dat.interp(x)))
HW.flow.wq.da$TP.kg=with(HW.flow.wq.da,Load.Calc.kg(value,TP.inter))
HW.flow.wq.da$TN.kg=with(HW.flow.wq.da,Load.Calc.kg(value,TN.inter))
HW.flow.wq.da=merge(HW.flow.wq.da,rbind(HW.flow.sites[,c("SITE","Region")],data.frame(SITE="S333R",Region="SRS")),"SITE")
HW.flow.wq.da$WY=WY(HW.flow.wq.da$Date.EST)
HW.flow.wq.da.region=ddply(HW.flow.wq.da,c("Region","Date.EST"),summarise,TFlow=sum(cfs.to.m3d(value),na.rm=T),TPLoad=sum(TP.kg,na.rm=T),TNLoad=sum(TN.kg,na.rm=T))
HW.flow.wq.da.region$WY=WY(HW.flow.wq.da.region$Date.EST)
quantile(subset(HW.flow.wq.da.region,Region=="SRS")$TFlow,probs=c(0.9))


## Hydrostats
HW.flow.wq.da.region$Q=HW.flow.wq.da.region$TFlow
HW.flow.wq.da.region$Date=HW.flow.wq.da.region$Date.EST
HW.flow.wq.da.region=merge(HW.flow.wq.da.region,analysis.periods,by.x="Date.EST",by.y="DATE",all.x=T)

range(subset(HW.flow.wq.da.region,Region=="SRS")$TFlow)
range(subset(HW.flow.wq.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL))$TFlow)
range(subset(HW.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL))$TFlow)

HighQ.all=quantile(subset(HW.flow.wq.da.region,Region=="SRS")$Q,probs=0.9,na.rm=T)
high.spells.v2(subset(HW.flow.wq.da.region,Region=="SRS")[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(HW.flow.wq.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL))[,c("Date","Q")],WY.type="FL",threshold=HighQ.all);
high.spells.v2(subset(HW.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL))[,c("Date","Q")],WY.type="FL",threshold=HighQ.all);

#ddply(HW.flow.wq.da.region,c("Region","WY"),summarise,HighQ=quantile(Q,probs=0.9,na.rm=T),HighDuration=mean(rle(as.numeric(Q>HighQ))$lengths[which(rle(as.numeric(Q>HighQ))$values==1)],na.rm=T))
#HighQ.all=quantile(subset(HW.flow.wq.da.region,Region=="SRS")$Q,probs=0.9,na.rm=T)
#ddply(HW.flow.wq.da.region,c("Region","WY"),summarise,HighDuration=mean(rle(as.numeric(Q>HighQ.all))$lengths[which(rle(as.numeric(Q>HighQ.all))$values==1)],na.rm=T))

range(subset(HW.flow.wq.da.region,Region=="TSPh")$TFlow)
range(subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(wilma.landfall.FL))$TFlow)
range(subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(irma.landfall.FL))$TFlow)

high.spells.v2(subset(HW.flow.wq.da.region,Region=="TSPh")[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(wilma.landfall.FL))[,c("Date","Q")],WY.type="FL");
high.spells.v2(subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(irma.landfall.FL))[,c("Date","Q")],WY.type="FL");

prepoststats=ddply(HW.flow.wq.da.region,c("Region","Storm","Period"),summarise,mean.val=format(mean(Q,na.rm=T),scientific = T),SD.val=format(sd(Q),scientific = T))
subset(prepoststats,Region=="TSPh"&Storm=="Wilma")
subset(prepoststats,Region=="TSPh"&Storm=="Irma")

boxplot(Q~Period+Storm+Region,subset(HW.flow.wq.da.region,is.na(Q)==F))
ddply(subset(HW.flow.wq.da.region,is.na(Storm)==F),c("Region","Storm"),summarise,chisq=kruskal.test(Q,as.factor(Period))$statistic,pval=kruskal.test(Q,as.factor(Period))$p.value,rounded.pval=round(pval,3))
test=cast(subset(HW.flow.wq.da.region,is.na(Storm)==F),Region+Storm~Period,value="Q",median)
test$eval=with(test,Pre<Post)
test

ddply(HW.flow.wq.da.region,c("Region"),summarise,val=RBFlash(Q))
Ann.RB=ddply(HW.flow.wq.da.region,c("Region","WY"),summarise,val=RBFlash(Q))
plot(val~WY,subset(Ann.RB,Region=="SRS"),type="l",ylim=c(0,1))
with(subset(Ann.RB,Region=="TSPh"),lines(WY,val,col="red"))
boxplot(val~Region,Ann.RB)
kruskal.test(val~Region,Ann.RB)

# WY Plot
HW.WY.load=ddply(HW.flow.wq.da,c("Region","WY"),summarise,TFlow=sum(cfs.to.m3d(value),na.rm=T),TPLoad=sum(TP.kg,na.rm=T),TNLoad=sum(TN.kg,na.rm=T))
HW.WY.load$TP.FWM=with(HW.WY.load,(TPLoad/TFlow)*1000000)
HW.WY.load$TN.FWM=with(HW.WY.load,(TNLoad/TFlow)*1000)

# TP Loads
ddply(HW.WY.load,c("Region"),summarise,min.val=min(TPLoad,na.rm=T),max.val=max(TPLoad,na.rm=T))
subset(HW.WY.load,Region=="SRS"&TPLoad==max(subset(HW.WY.load,Region=="SRS")$TPLoad))
subset(HW.WY.load,Region=="SRS"&TPLoad==max(subset(HW.WY.load,Region=="SRS"&WY!=2006)$TPLoad))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(wilma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS")$TP.kg))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(irma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS")$TP.kg))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(wilma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(irma.landfall.FL))$TP.kg))

subset(HW.WY.load,Region=="TSPh"&TPLoad==max(subset(HW.WY.load,Region=="TSPh")$TPLoad))
subset(HW.WY.load,Region=="TSPh"&TPLoad==max(subset(HW.WY.load,Region=="TSPh"&WY!=2018)$TPLoad))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(wilma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh")$TP.kg))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(irma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh")$TP.kg))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(wilma.landfall.FL))$TP.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(irma.landfall.FL))$TP.kg))

# TN Loads
ddply(HW.WY.load,c("Region"),summarise,min.val=min(TNLoad,na.rm=T),max.val=max(TNLoad,na.rm=T))
subset(HW.WY.load,Region=="SRS"&TNLoad==max(subset(HW.WY.load,Region=="SRS")$TNLoad))
subset(HW.WY.load,Region=="SRS"&TNLoad==max(subset(HW.WY.load,Region=="SRS"&WY!=2006)$TNLoad))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(wilma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS")$TN.kg))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(irma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS")$TN.kg))
ks.test(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(wilma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="SRS"&WY==WY(irma.landfall.FL))$TN.kg))

subset(HW.WY.load,Region=="TSPh"&TNLoad==max(subset(HW.WY.load,Region=="TSPh")$TNLoad))
subset(HW.WY.load,Region=="TSPh"&TNLoad==max(subset(HW.WY.load,Region=="TSPh"&WY!=2006)$TNLoad))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(wilma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh")$TN.kg))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(irma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh")$TN.kg))
ks.test(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(wilma.landfall.FL))$TN.kg,ecdf(subset(HW.flow.wq.da,Region=="TSPh"&WY==WY(irma.landfall.FL))$TN.kg))


#Supplemental Figure
regions.val=c("SRS","TSPh")
regions.labels=c("SRS","TS/Ph")
ylim.val=c(0,1.05);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"FigS2_LoadCDF.tiff"),width=6.5,height=5.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(3,2,0.5,1),oma=c(2,3,1,0.5));
layout(matrix(c(1:4,5,5),3,2,byrow=T),heights=c(rep(1,2),0.25))

xlim.val.max=c(46,11)
by.x.val=c(10,2)
for(i in 1:2){
  xlim.val=c(0,xlim.val.max[i]);by.x=by.x.val[i];xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
  plot(proportion~value,ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i])$TP.kg),type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,yaxs="i",xaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf.CI(subset(HW.flow.wq.da,Region==regions.val[i])$TP.kg),shaded.range(value,lwr,upr,"grey",lty=1))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i])$TP.kg),lines(value,proportion,col="grey",lwd=1.5))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i]&WY==WY(wilma.landfall.FL))$TP.kg),lines(value,proportion,col="indianred1",lwd=1.5,lty=2))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i]&WY==WY(irma.landfall.FL))$TP.kg),lines(value,proportion,col="dodgerblue1",lwd=1.5))
  axis_fun(1,line=-0.25,xmaj,xmin,xmaj,1);axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
  mtext(side=3,regions.labels[i])
  if(i==1){mtext(side=2,line=2.5,"Proportion")}
  mtext(side=1,line=2,"TP Load (kg d\u207B\u00B9)")
}
xlim.val.max=c(5500,1250)
by.x.val=c(1e3,0.25e3)
for(i in 1:2){
  xlim.val=c(0,xlim.val.max[i]);by.x=by.x.val[i];xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
  plot(proportion~value,ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i])$TN.kg),type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,yaxs="i",xaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf.CI(subset(HW.flow.wq.da,Region==regions.val[i])$TN.kg),shaded.range(value,lwr,upr,"grey",lty=1))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i])$TN.kg),lines(value,proportion,col="grey",lwd=1.5))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i]&WY==WY(wilma.landfall.FL))$TN.kg),lines(value,proportion,col="indianred1",lwd=1.5,lty=2))
  with(ecdf.v2(subset(HW.flow.wq.da,Region==regions.val[i]&WY==WY(irma.landfall.FL))$TN.kg),lines(value,proportion,col="dodgerblue1",lwd=1.5))
  axis_fun(1,line=-0.5,xmaj,xmin,xmaj,1);axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
  mtext(side=1,line=2,"TN Load (kg d\u207B\u00B9)")
  if(i==1){mtext(side=2,line=2.5,"Proportion")}
}
plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
legend.text=c("POR (\u00B1 95% Confidence Interval)",paste0("WY",WY(wilma.landfall.FL)," (Wilma)"),paste0("WY",WY(irma.landfall.FL)," (Irma)"))
pt.cols=c(adjustcolor("grey",0.25),"indianred1","dodgerblue1")
legend(0.5,0.8,legend=legend.text,lty=c(1,2,1),col=c("grey","indianred1","dodgerblue1"),lwd=c(1,1.5,1.5),pch =c(22,NA,NA),pt.bg=pt.cols,pt.cex=2,ncol=1,cex=1.25,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()


plot(TFlow~WY,HW.WY.load,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
with(subset(HW.WY.load,Region=="TSPh"),shaded.range(WY,rep(-10,length(WY)),TFlow,"grey",lty=1))
with(cast(HW.WY.load,WY~Region,value="TFlow",mean),shaded.range(WY,TSPh,TSPh+Tide,"grey50",lty=1))

regions.val=c("SRS","TSPh")
regions.labels=c("SRS","TS/Ph")
xlim.val=c(2005,2018);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x);
ylab.line=2.75
#tiff(filename=paste0(plot.path,"Fig6_InflowFlowLoad.tiff"),width=6.5,height=5.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3.25,0.5,1.5),oma=c(2,3.25,1,3.25),mgp=c(3,1,0));
layout(matrix(c(1:4,5:7,4),4,2,byrow=F),heights=c(rep(1,3),0.25))

for(i in 1:2){
  ylim.val=c(0,2.1e9);by.y=5e8;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(HW.WY.load,Region==regions.val[i]),shaded.range(WY,rep(-10,length(WY)),TFlow,"grey",lty=1))
  if(i==2){with(cast(HW.WY.load,WY~Region,value="TFlow",mean),shaded.range(WY,TSPh,TSPh+Tide,"grey50",lty=1))}
  text(x=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)-0.25),y=rep(ylim.val[2]-1e8,2),label=c("Wilma","Irma"))
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  axis_fun(2,ymaj,ymin,ymaj/1e8,1);axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"Discharge\n(x10\u2078; m\u00B3 Yr\u207B\u00B9)")}
  mtext(side=3,regions.labels[i])
  
  ylim.val=c(0,15e3);by.y=5e3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(HW.WY.load,Region==regions.val[i]),shaded.range(WY,rep(-10,length(WY)),TPLoad,"indianred1",lty=1))
  #if(i==2){with(cast(HW.WY.load,WY~Region,value="TPLoad",mean),shaded.range(WY,TSPh,TSPh+Tide,"indianred4",lty=1))}
  axis_fun(2,ymaj,ymin,format(ymaj/1e3),1);axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"TP Load\n(x10\u00B3; kg Yr\u207B\u00B9)")}
  ylim.val=c(0,15);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  with(subset(HW.WY.load,Region==regions.val[i]),pt_line(WY,TP.FWM,1,"black",1,21,"indianred1"))
  axis_fun(4,ymaj,ymin,ymaj,1)
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  if(i==2){mtext(side=4,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",las=0)}
  
  ylim.val=c(0,200e4);by.y=500e3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(HW.WY.load,Region==regions.val[i]),shaded.range(WY,rep(-10,length(WY)),TNLoad,"dodgerblue",lty=1))
  #if(i==2){with(cast(HW.WY.load,WY~Region,value="TNLoad",mean),shaded.range(WY,TSPh,TSPh+Tide,"dodgerblue4",lty=1))}
  axis_fun(2,ymaj,ymin,format(ymaj/1e3),1);axis_fun(1,line=-0.5,xmaj,xmin,xmaj,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"TN Load\n(x10\u00B3; kg Yr\u207B\u00B9)")}
  ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  with(subset(HW.WY.load,Region==regions.val[i]),pt_line(WY,TN.FWM,1,"black",1,21,"dodgerblue"))
  axis_fun(4,ymaj,ymin,format(ymaj),1)
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  if(i==2){mtext(side=4,line=2,"TN FWM (mg L\u207B\u00B9)",las=0)}
  mtext(side=1,line=2,"Water Year")
  
  if(i==1){
    plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
    legend.text=c("Flow","TP Load","TP FWM","Flow to Tide ","TN Load","TN FWM")
    pt.cols=c(adjustcolor("grey",0.25),adjustcolor("indianred1",0.25),"indianred1",adjustcolor("grey50",0.25),adjustcolor("dodgerblue1",0.25),"dodgerblue1")
    legend(0.5,0.8,legend=legend.text,lty=c(0,0,1,0,0,1),col=c("grey","indianred1","black","grey50","dodgerblue1","black"),lwd=c(1,1,1,1,1,1),pch =c(22,22,21,22,22,21),pt.bg=pt.cols,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
  }
}
dev.off()

###
estuary.flow.sites=data.frame(Station=c("02290829","BROADUPS"),
                              RiverName=c("Bottle","BroadUp"),
                              DBKEY=c("89858","90028"),
                              WQSite=c("P36","P36"))

estuary.wq=SFWMD.DBHYDRO.Data.WQ(sdate,edate,"P36",wq.params$Test.Number)
estuary.wq=subset(estuary.wq,Collection.Method=="G")
estuary.wq$AbsMDL=with(estuary.wq,ifelse(Value<0,abs(Value),Value))
estuary.wq=merge(estuary.wq,wq.params,"Test.Number")
estuary.wq.xtab=cast(estuary.wq,Station.ID+Date.EST~param,value="AbsMDL",mean)
estuary.wq.xtab$TN.final=with(estuary.wq.xtab,SFWMD.TN.Combine(NOx,TKN,TN))

estuary.q=data.frame()
for(i in 1:nrow(estuary.flow.sites)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,estuary.flow.sites$DBKEY[i])
  tmp$DBKEY=as.character(estuary.flow.sites$DBKEY[i])
  tmp$Station=as.character(estuary.flow.sites$Station[i])
  estuary.q=rbind(estuary.q,tmp)
  print(i)
}
estuary.q=merge(estuary.q,estuary.flow.sites,c("Station","DBKEY"))
estuary.q$Date.EST=date.fun(estuary.q$Date)
estuary.q=merge(estuary.q,analysis.periods,by.x="Date.EST",by.y="DATE",all.x=T)
estuary.q$WY=WY(estuary.q$Date.EST)
estuary.q$pos.flow.cfs=with(estuary.q,ifelse(Data.Value<0,0,Data.Value))
estuary.q$neg.flow.cfs=with(estuary.q,ifelse(Data.Value<0,abs(Data.Value),0))
estuary.q$Q.m3d=cfs.to.m3d(estuary.q$Data.Value)
estuary.q$Q.m3d.pos=with(estuary.q,ifelse(Q.m3d<0,0,Q.m3d))


plot(Q.m3d.pos~Date.EST,subset(estuary.q,RiverName=="BroadUp"))

vars=c("Date.EST","WY","RiverName","WQSite","pos.flow.cfs","neg.flow.cfs")
estuary.flow.wq.da=merge(estuary.q[,vars],estuary.wq.xtab,by.x=c("Date.EST","WQSite"),by.y=c("Date.EST","Station.ID"),all.x=T)
estuary.flow.wq.da=estuary.flow.wq.da[order(estuary.flow.wq.da$RiverName,estuary.flow.wq.da$Date.EST),]
estuary.flow.wq.da$TP.int=with(estuary.flow.wq.da,ave(TP,RiverName,FUN=function(x)dat.interp(x)))
estuary.flow.wq.da$TN.int=with(estuary.flow.wq.da,ave(TN.final,RiverName,FUN=function(x)dat.interp(x)))
estuary.flow.wq.da$TP.kg=with(estuary.flow.wq.da,Load.Calc.kg(pos.flow.cfs,TP.int))
estuary.flow.wq.da$TN.kg=with(estuary.flow.wq.da,Load.Calc.kg(pos.flow.cfs,TN.int))

estuary.wy=ddply(estuary.flow.wq.da,c("RiverName","WY"),summarise,TFlow=sum(cfs.to.m3d(pos.flow.cfs),na.rm=T),TPLoad=sum(TP.kg,na.rm=T),TNLoad=sum(TN.kg,na.rm=T),N.flow=N(pos.flow.cfs))
estuary.wy$TP.FWM=with(estuary.wy,(TPLoad/TFlow)*1000000)
estuary.wy$TN.FWM=with(estuary.wy,(TNLoad/TFlow)*1000)

### 
coastal.flow.sites=data.frame(Station=c("HARNEYFLAM","SHARKRIVBG","02290878","TAYLORS3","EASTCK","MUD_CRKM"),
                              RiverName=c("Harney","Shark","Broad","Taylor","East","Mud"),
                              Region=c(rep("SRS",3),rep("TS",3)),
                              DBKEY=c("90144","90243","89871","AN865","AM250","AN862"),
                              WQSite=c("FLAB37","FLAB39","FLAB33",rep("TS/PH7a",3)))

coastal.wq.sites=data.frame(RiverName=c("Harney","Shark","Broad","Taylor","East","Mud"),
                            WQSite=c("FLAB37","FLAB39","FLAB33",rep("TS/PH7a",3)),
                            source=c(rep("SFWMD",3),rep("FCE",3)))
coast.wq=SFWMD.DBHYDRO.Data.WQ(sdate,edate,subset(coastal.wq.sites,source=="SFWMD")$WQSite,wq.params$Test.Number)
coast.wq=subset(coast.wq,Collection.Method=="G")
coast.wq$AbsMDL=with(coast.wq,ifelse(Value<0,abs(Value),Value))
coast.wq=merge(coast.wq,wq.params,"Test.Number")
coast.wq=merge(coast.wq,coastal.wq.sites,by.x="Station.ID",by.y="WQSite")
coast.wq.xtab=cast(coast.wq,Station.ID+RiverName+Date.EST~param,value="AbsMDL",mean)
coast.wq.xtab$TN.final=with(coast.wq.xtab,SFWMD.TN.Combine(NOx,TKN,TN))

plot(TP~Date.EST,subset(coast.wq.xtab,Station.ID=="FLAB39"),type="b")
mtext(side=2,line=2,"TP");mtext(side=3,"FLAB37")

coast.wq.xtab2=read.csv(paste0(data.path,"SFWMD/CWQMN_thru_April2018_For_Analysis.csv"))
coast.wq.xtab2$Date.EST=date.fun(coast.wq.xtab2$Collection.Date,form="%m/%d/%Y")
coast.wq.xtab2$TP.mgL=(coast.wq.xtab2$CTP.P.um.l*P.mw)*0.001
coast.wq.xtab2$TN.mgL=coast.wq.xtab2$TN.mg.l
coast.wq.xtab2$TOC.mgL=(coast.wq.xtab2$TOC.um.l*C.mw)*0.001
coast.wq.xtab2$DOC.mgL=NA
coast.wq.xtab2=coast.wq.xtab2[,c("Date.EST","Station.ID","TP.mgL","TN.mgL","TOC.mgL","DOC.mgL")]

plot(TOC.mgL~Date.EST,subset(coast.wq.xtab2,Station.ID=="FLAB39"),type="b")
plot(TN.mgL~Date.EST,subset(coast.wq.xtab2,Station.ID=="FLAB37"),type="b")
plot(TP.mgL~Date.EST,subset(coast.wq.xtab2,Station.ID=="FLAB37"),type="b")

losada=read.csv(paste0(data.path,"/WQ/LT_ND_Losada_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
losada$Date.EST=date.fun(as.character(losada$DATE),form="%Y-%m-%d")
losada$TP.mgL=(losada$TP*P.mw)*0.001
losada$TN.mgL=(losada$TN*N.mw)*0.001
losada$TOC.mgL=(losada$TOC*C.mw)*0.001
losada$DOC.mgL=(losada$DOC*C.mw)*0.001
losada$Station.ID=losada$SITENAME
losada2=losada[,c("Date.EST","Station.ID","TP.mgL","TN.mgL","TOC.mgL","DOC.mgL")]

coast.wq.xtab2=subset(rbind(coast.wq.xtab2,losada2),Station.ID%in%coastal.wq.sites$WQSite)  
coast.wq.xtab2$TOC.mgL2=with(coast.wq.xtab2,ifelse(is.na(TOC.mgL)==T,DOC.mgL,TOC.mgL))

plot(TOC.mgL2~Date.EST,subset(coast.wq.xtab2,Station.ID=="TS/PH7a"))
plot(TOC.mgL2~Date.EST,subset(coast.wq.xtab2,Station.ID=="FLAB37"))

plot(DOC~TOC,losada);abline(0,1)
lm(DOC~TOC,losada)
## Use FCE grab data for TS load TS/Ph 7a or 7b
coastal.q=data.frame()
for(i in 1:nrow(coastal.flow.sites)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,coastal.flow.sites$DBKEY[i])
  tmp$DBKEY=as.character(coastal.flow.sites$DBKEY[i])
  tmp$Station=as.character(coastal.flow.sites$Station[i])
  coastal.q=rbind(coastal.q,tmp)
  print(i)
}
coastal.q=merge(coastal.q,coastal.flow.sites,c("Station","DBKEY"))
coastal.q$Date.EST=date.fun(coastal.q$Date)
coastal.q=merge(coastal.q,analysis.periods,by.x="Date.EST",by.y="DATE",all.x=T)
coastal.q$WY=WY(coastal.q$Date.EST)
coastal.q$hydroday=hydro.day(coastal.q$Date.EST,"FL")
coastal.q$pos.flow.cfs=with(coastal.q,ifelse(Data.Value<0,0,Data.Value))
coastal.q$neg.flow.cfs=with(coastal.q,ifelse(Data.Value<0,abs(Data.Value),0))
coastal.q$Q.m3d=cfs.to.m3d(coastal.q$Data.Value)
coastal.q$Q.m3d.pos=with(coastal.q,ifelse(Q.m3d<0,0,Q.m3d))
coastal.q$month=format(coastal.q$Date.EST,"%m")
coastal.q$CY=format(coastal.q$Date.EST,"%Y")

vars=c("Date.EST","WY","RiverName","WQSite","Region","pos.flow.cfs","neg.flow.cfs")
coast.flow.wq.da=merge(coastal.q[,vars],coast.wq.xtab2,by.x=c("Date.EST","WQSite"),by.y=c("Date.EST","Station.ID"),all.x=T)
coast.flow.wq.da=coast.flow.wq.da[order(coast.flow.wq.da$RiverName,coast.flow.wq.da$Date.EST),]
coast.flow.wq.da$TP.int=with(coast.flow.wq.da,ave(TP.mgL,RiverName,FUN=function(x)dat.interp(x)))
coast.flow.wq.da$TN.int=with(coast.flow.wq.da,ave(TN.mgL,RiverName,FUN=function(x)dat.interp(x)))
coast.flow.wq.da$TOC.int=with(coast.flow.wq.da,ave(TOC.mgL2,RiverName,FUN=function(x)dat.interp(x)))
coast.flow.wq.da$TP.kg=with(coast.flow.wq.da,Load.Calc.kg(pos.flow.cfs,TP.int))
coast.flow.wq.da$TN.kg=with(coast.flow.wq.da,Load.Calc.kg(pos.flow.cfs,TN.int))
coast.flow.wq.da$TOC.kg=with(coast.flow.wq.da,Load.Calc.kg(pos.flow.cfs,TOC.int))

coast.flow.wq.da$TP.kg.landward=with(coast.flow.wq.da,Load.Calc.kg(neg.flow.cfs,TP.int))
coast.flow.wq.da$TN.kg.landward=with(coast.flow.wq.da,Load.Calc.kg(neg.flow.cfs,TN.int))
coast.flow.wq.da$TOC.kg.landward=with(coast.flow.wq.da,Load.Calc.kg(neg.flow.cfs,TOC.int))

coast.load.da.region=ddply(coast.flow.wq.da,c("Region","Date.EST","WY"),summarise,GOM.TP.kg=sum(TP.kg,na.rm=T),ENP.TP.kg=sum(TP.kg.landward,na.rm=T))
plot(GOM.TP.kg~Date.EST,subset(coast.load.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL)),type="l")
with(subset(coast.load.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL)),lines(Date.EST,ENP.TP.kg,col="Red"));
abline(v=wilma.landfall.FL)
plot(GOM.TP.kg~Date.EST,subset(coast.load.da.region,Region=="SRS"&WY==WY(irma.landfall.FL)),type="l")
with(subset(coast.load.da.region,Region=="SRS"&WY==WY(irma.landfall.FL)),lines(Date.EST,ENP.TP.kg,col="Red"));
abline(v=irma.landfall.FL)

## Hydrostats
coast.flow.wq.da$Q=coast.flow.wq.da$pos.flow.cfs 
coast.flow.wq.da$Date=coast.flow.wq.da$Date.EST
coast.flow.wq.da=merge(coast.flow.wq.da,analysis.periods,by.x="Date.EST",by.y="DATE",all.x=T)

plot(Q~Date,subset(coast.flow.wq.da,RiverName=="Shark"),pch=21,bg="grey")
with(subset(coast.flow.wq.da,RiverName=="Harney"),points(Date,Q,pch=21,bg=adjustcolor("indianred1",0.25)))
with(subset(coast.flow.wq.da,RiverName=="Broad"),points(Date,Q,pch=21,bg=adjustcolor("dodgerblue1",0.25)))
tmp=subset(coast.flow.wq.da,Region=="SRS")
tmp$RiverName=factor(tmp$RiverName,levels=c("Broad","Harney","Shark"))
boxplot(Q~RiverName,tmp,outline=F);kruskal.test(Q~RiverName,tmp)

coast.flow.wq.da$RiverName=factor(coast.flow.wq.da$RiverName,levels=c("Broad","Harney","Shark","Taylor","East","Mud"))

#tiff(filename=paste0(plot.path,"FigS1_dailydischarge.tiff"),width=6,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,1.5,0.5,1.5),oma=c(2,2,0.5,0.5),mgp=c(3,1,0));

ylim.val=c(0,4.1e6);by.y=1e6;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
boxplot(cfs.to.m3d(Q)~RiverName,coast.flow.wq.da,outline=F,ylim=ylim.val,yaxt="n",xaxt="n",col="grey")
abline(v=3.5,lty=1)
text(c(2,5),rep(ylim.val[2]+0.05e6,2),c("SRS","TS/Ph"),cex=0.9)
axis_fun(2,ymaj,ymin,ymaj/1e5,1)
axis_fun(1,line=-0.5,1:6,1:6,c("Broad","Harney","Shark","Taylor","East","Mud"),1)
mtext(side=2,line=2,"Daily Discharge (x10\u2075 m\u00B3 d\u207B\u00B9)")
mtext(side=1,line=2,"Coastal River Discharge")
dev.off()

coast.flow.wq.da.region=ddply(coast.flow.wq.da,c("Region","Date","WY","Storm","Period"),summarise,Q=sum(cfs.to.m3d(pos.flow.cfs),na.rm=T))

#SRS summary statsitics
format(range(subset(coast.flow.wq.da.region,Region=="SRS")$Q),scientific=T)
range(subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL))$Q)
range(subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL))$Q)

plot(Q~Date,subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL)));abline(v=wilma.landfall.FL)
plot(Q~Date,subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL)));abline(v=irma.landfall.FL)

high.spells.v2(subset(coast.flow.wq.da.region,Region=="SRS")[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(wilma.landfall.FL))[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(coast.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL))[,c("Date","Q")],WY.type="FL")

#TS summary statsitics
cast(coast.flow.wq.da,WY~RiverName,value="Q",fun.aggregate=function(x)sum(x,na.rm=T))
TS.coast=subset(coast.flow.wq.da.region,Region=="TS"&WY%in%seq(2007,2018,1))

format(range(TS.coast$Q),scientific=T)
subset(TS.coast,Q==max(TS.coast$Q,na.rm=T))
plot(Q~Date,subset(coast.flow.wq.da.region,Region=="TS"&WY==2012),type="l")
range(subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(wilma.landfall.FL))$Q)
range(subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(irma.landfall.FL))$Q)

plot(Q~Date,subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(wilma.landfall.FL)));abline(v=wilma.landfall.FL)
plot(Q~Date,subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(irma.landfall.FL)));abline(v=irma.landfall.FL)

high.spells.v2(TS.coast[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(wilma.landfall.FL))[,c("Date","Q")],WY.type="FL")
high.spells.v2(subset(coast.flow.wq.da.region,Region=="TS"&WY==WY(irma.landfall.FL))[,c("Date","Q")],WY.type="FL")



ddply(subset(coast.flow.wq.da.region,is.na(Storm)==F),c("Region","Storm"),summarise,chisq=kruskal.test(Q,as.factor(Period))$statistic,pval=kruskal.test(Q,as.factor(Period))$p.value,rounded.pval=round(pval,3))
prepoststats=ddply(coast.flow.wq.da.region,c("Region","Storm","Period"),summarise,mean.val=mean(Q,na.rm=T),SD.val=sd(Q,na.rm=T))
subset(prepoststats,Region=="SRS"&Storm=="Wilma")
subset(prepoststats,Region=="SRS"&Storm=="Irma")
subset(prepoststats,Region=="TS"&Storm=="Wilma")
subset(prepoststats,Region=="TS"&Storm=="Irma")

# WY Plot
coastal.wy=ddply(coast.flow.wq.da,c("RiverName","Region","WY"),summarise,TFlow=sum(cfs.to.m3d(pos.flow.cfs),na.rm=T),TPLoad=sum(TP.kg,na.rm=T),TNLoad=sum(TN.kg,na.rm=T),TOCLoad=sum(TOC.kg,na.rm=T),N.flow=N(pos.flow.cfs))
coastal.wy$TP.FWM=with(coastal.wy,(TPLoad/TFlow)*1000000)
coastal.wy$TN.FWM=with(coastal.wy,(TNLoad/TFlow)*1000)
coastal.wy$TOC.FWM=with(coastal.wy,(TOCLoad/TFlow)*1000)

## Trend analysis
coastal.wy.region=ddply(coast.flow.wq.da,c("Region","WY"),summarise,TFlow=sum(cfs.to.m3d(pos.flow.cfs),na.rm=T),TPLoad=sum(TP.kg,na.rm=T),TNLoad=sum(TN.kg,na.rm=T),TOCLoad=sum(TOC.kg,na.rm=T),N.flow=N(pos.flow.cfs))
coastal.wy.region$TP.FWM=with(coastal.wy.region,(TPLoad/TFlow)*1000000)
coastal.wy.region$TN.FWM=with(coastal.wy.region,(TNLoad/TFlow)*1000)
coastal.wy.region$TOC.FWM=with(coastal.wy.region,(TOCLoad/TFlow)*1000)

##Quick compare of FW and Coastal TP Loads
test=merge(merge(HW.WY.load,data.frame(Region=c("SRS","Tide","TSPh"),Region2=c("SRS",NA,"TS")),"Region"),
           coastal.wy.region,by.x=c("Region2","WY"),by.y=c("Region","WY"))

plot(TNLoad.y~WY,subset(test,Region=="SRS"),ylim=c(10,4000000),log="y");#coastal
with(subset(test,Region=="SRS"),lines(WY,TNLoad.x));#FW
with(subset(estuary.wy,RiverName=="Bottle"),lines(WY,TNLoad,col="red"))
with(subset(estuary.wy,RiverName=="BroadUp"),lines(WY,TNLoad,col="red",lty=2))
with(ddply(estuary.wy,"WY",summarise,load=sum(TNLoad,na.rm=T)),lines(WY,load,col="blue"))
###



#srs.coast.load=subset(coastal.wy.region,Region=="SRS")
#srs.coast.load$Cum.TPLoad=cumsum(srs.coast.load$TPLoad)
#srs.coast.load$Cum.TNLoad=cumsum(srs.coast.load$TNLoad)
#srs.coast.load$Cum.TOCLoad=cumsum(srs.coast.load$TOCLoad)
#plot(Cum.TPLoad~WY,srs.coast.load)
#plot(Cum.TNLoad~WY,srs.coast.load)
#plot(Cum.TOCLoad~WY,srs.coast.load)
####
#test=lm(TOCLoad~WY,subset(coastal.wy.region,Region=="SRS"))
#test.seg=segmented(test,seg.Z=~WY)
#summary(test.seg)
#plot(TOCLoad~WY,subset(coastal.wy.region,Region=="SRS"))
#plot(test.seg,add=T)
####

#Compare cumulative distributions
# TP Loads
ddply(coastal.wy.region,c("Region"),summarise,min.val=min(TPLoad,na.rm=T),max.val=max(TPLoad,na.rm=T))
range(subset(coastal.wy.region,Region=="SRS")$TPLoad)
subset(coastal.wy.region,Region=="SRS"&TPLoad==max(subset(coastal.wy.region,Region=="SRS")$TPLoad))

range(subset(coastal.wy.region,Region=="TS"&WY>2006)$TPLoad)
subset(coastal.wy.region,Region=="TS"&WY>2006&TPLoad==max(subset(coastal.wy.region,Region=="TS"&WY>2006)$TPLoad))
test=ddply(subset(coastal.wy,RiverName%in%c("Taylor","Mud")),"WY",summarise,load=sum(TPLoad,na.rm=T))
range(test$load)
subset(test,load==max(test$load))

ddply(coastal.wy.region,c("Region"),summarise,min.val=min(TNLoad,na.rm=T),max.val=max(TNLoad,na.rm=T))
range(subset(coastal.wy.region,Region=="SRS")$TNLoad)
subset(coastal.wy.region,Region=="SRS"&TPLoad==max(subset(coastal.wy.region,Region=="SRS")$TPLoad))


coastal.wy$RiverName=factor(coastal.wy$RiverName,levels=c("Broad River","Harney","Shark","Taylor","Mud","East"))
regions.val=c("SRS","TS")
regions.labels=c("SRS","TS/Ph")
ylab.line=2.75
xlim.val=c(2005,2018);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x);
#tiff(filename=paste0(plot.path,"Fig7_CoastalFlowLoad.tiff"),width=6.5,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,3,0.5,1.5),oma=c(1.5,3.25,1,3.25),mgp=c(3,1,0));
layout(matrix(c(1:10),5,2,byrow=F),heights=c(rep(1,4),0.3))

for(i in 1:2){
  ylim.val=list(c(0,1.6e9),c(0,2e8))[[i]];by.y=list(0.5e9,0.5e8)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TFlow~WY,coastal.wy,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  flow.tmp=cast(subset(coastal.wy,Region==regions.val[i]),WY~RiverName,value="TFlow",mean)
  flow.tmp[is.na(flow.tmp)]=0
  stacked.plot.poly(flow.tmp[,1],data.matrix(flow.tmp[,-1]),c("grey80","grey50","grey30"),lty=1)
  text(x=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)-0.25),y=list(ylim.val[2]-0.1e9,ylim.val[2]-0.1e8)[[i]],label=c("Wilma","Irma"))
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  axis_fun(2,ymaj,ymin,format(ymaj/1e8),1);axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"Discharge\n(x10\u2079; m\u00B3 Yr\u207B\u00B9)")}
  mtext(side=3,regions.labels[i])
  
  ylim.val=list(c(0,40e3),c(0,5e3))[[i]];by.y=list(10e3,1e3)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TPLoad~WY,coastal.wy,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  TPload.tmp=cast(subset(coastal.wy,Region==regions.val[i]),WY~RiverName,value="TPLoad",mean)
  TPload.tmp[is.na(TPload.tmp)]=0
  stacked.plot.poly(TPload.tmp[,1],data.matrix(TPload.tmp[,-1]),c("grey80","grey50","grey30"),lty=1)
  axis_fun(2,ymaj,ymin,format(ymaj/1e3),1);axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"TP Load\n(x10\u00B3; kg Yr\u207B\u00B9)")}
  fwm=ddply(subset(coastal.wy,Region==regions.val[i]),c("Region","WY"),summarise,Flow=sum(TFlow),Load=sum(TPLoad),FWM=(Load/Flow)*1000000)
  ylim.val=list(c(0,40),c(0,30))[[i]];by.y=list(10,10)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  with(fwm,pt_line(WY,FWM,1,"grey50",1,21,"indianred1",cex=1.25,pt.lwd = 0.1))
  axis_fun(4,ymaj,ymin,ymaj,1)
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  if(i==2){mtext(side=4,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",las=0)}
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  
  ylim.val=list(c(0,1250e3),c(0,150e3))[[i]];by.y=list(500e3,50e3)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TPLoad~WY,coastal.wy,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  TNload.tmp=cast(subset(coastal.wy,Region==regions.val[i]),WY~RiverName,value="TNLoad",mean)
  TNload.tmp[is.na(TNload.tmp)]=0
  stacked.plot.poly(TNload.tmp[,1],data.matrix(TNload.tmp[,-1]),c("grey80","grey50","grey30"),lty=1)
  axis_fun(2,ymaj,ymin,format(ymaj/1e3),1);axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"TN Load\n(x10\u00B3; kg Yr\u207B\u00B9)")}
  fwm=ddply(subset(coastal.wy,Region==regions.val[i]),c("Region","WY"),summarise,Flow=sum(TFlow),Load=sum(TNLoad),FWM=(Load/Flow)*1000)
  ylim.val=c(0,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  with(fwm,pt_line(WY,FWM,1,"grey50",1,21,"indianred1",cex=1.25,pt.lwd = 0.1))
  axis_fun(4,ymaj,ymin,format(ymaj),1)
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  if(i==2){mtext(side=4,line=2,"TN FWM (mg L\u207B\u00B9)",las=0)}
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  
  ylim.val=list(c(0,21e6),c(0,2.5e6))[[i]];by.y=list(5e6,0.5e6)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TPLoad~WY,coastal.wy,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  TOCload.tmp=cast(subset(coastal.wy,Region==regions.val[i]),WY~RiverName,value="TOCLoad",mean)
  TOCload.tmp[is.na(TOCload.tmp)]=0
  stacked.plot.poly(TOCload.tmp[,1],data.matrix(TOCload.tmp[,-1]),c("grey80","grey50","grey30"),lty=1)
  axis_fun(2,ymaj,ymin,format(ymaj/1e6),1);axis_fun(1,line=-0.5,xmaj,xmin,xmaj,1);box(lwd=1)
  if(i==1){mtext(side=2,line=ylab.line,"TOC Load\n(x10\u2076; kg Yr\u207B\u00B9)")}
  fwm=ddply(subset(coastal.wy,Region==regions.val[i]),c("Region","WY"),summarise,Flow=sum(TFlow),Load=sum(TOCLoad),FWM=(Load/Flow)*1000)
  ylim.val=c(0,20);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(TFlow~WY,HW.WY.load,yaxt="n",xaxt="n",ylab=NA,xlab=NA,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i")
  with(fwm,pt_line(WY,FWM,1,"grey50",1,21,"indianred1",cex=1.25,pt.lwd = 0.1))
  axis_fun(4,ymaj,ymin,ymaj,1)
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  if(i==2){mtext(side=4,line=2,"TOC FWM (mg L\u207B\u00B9)",las=0)}
  abline(v=c(WY(wilma.landfall.FL),WY(irma.landfall.FL)),col="black",lty=2)
  mtext(side=1,line=1.75,"Water Year")
  
  plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
  legend.text=list(c("Broad","Harney","Shark","FWM"),c("Taylor","Mud","East","FWM"))[[i]]
  pt.cols=c(adjustcolor(c("grey80","grey50","grey30"),0.25),"indianred1")
  legend(0.5,0.2,legend=legend.text,lty=c(0,0,0,1),col=c(c("grey80","grey50","grey30"),"grey50"),lwd=c(1,1,1,1,1,1),pch=c(22,22,22,21),pt.bg=pt.cols,pt.cex=1.5,ncol=4,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
}
dev.off()






