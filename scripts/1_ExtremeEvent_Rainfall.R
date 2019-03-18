##

# Rainfall ----------------------------------------------------------------
rf.mon=subset(sfwmd.mon,ACTIVITY_S=='Rain')

enp.buffer=gBuffer(enp,width=5*1000)
enp.buffer=SpatialPolygonsDataFrame(enp.buffer,data.frame(row.names = "buffer",width.km=5))

rf.mon.enp=rf.mon[enp.buffer,]
paste(rf.mon.enp@data$STATION,collapse="/")

dbkey.list=read.csv(paste0(data.path,"RF_DBKeys_REPORT.csv"))
dbkey.list$Start.Date2=date.fun(as.POSIXct(as.character(dbkey.list$Start.Date),format="%d-%B-%Y"))
dbkey.list$End.Date2=date.fun(as.POSIXct(as.character(dbkey.list$End.Date),format="%d-%B-%Y"))
dbkey.list$diff.Date=with(dbkey.list,as.numeric((End.Date2-Start.Date2)/365))
dbkey.list$Start.Year=as.numeric(format(dbkey.list$Start.Date2,"%Y"))
dbkey.list$End.Year=as.numeric(format(dbkey.list$End.Date2,"%Y"))

nchar(as.character(dbkey.list$Dbkey))
range(dbkey.list$Start.Year,na.rm=T)
dbkey.list=subset(dbkey.list,Start.Year<=1999&End.Year>=2019)
#dbkey.list=subset(dbkey.list,diff.Date>20&Start.Year>1992)
ddply(dbkey.list,c("Station"),summarise,N.val=N(Dbkey))
plot(subset(rf.mon.enp,STATION%in%ddply(dbkey.list,c("Station"),summarise,N.val=N(Dbkey))$Station));plot(enp,add=T)

rf.mon.enp=subset(rf.mon.enp,STATION%in%ddply(dbkey.list,c("Station"),summarise,N.val=N(Dbkey))$Station)
rf.mon.enp.th=thessian_create.v2(rf.mon.enp,enp,plot=F)
rf.mon.enp.th$area.pro=with(rf.mon.enp.th@data,area/sum(area,na.rm=T))

tm_shape(rf.mon.enp.th)+tm_polygons(col="dodgerblue2",alpha=0.2)+
  tm_shape(rf.mon.enp)+tm_dots()

#Download rainfall
sdate=as.Date("1999-05-01");# Begining of Florida WY2000
edate=as.Date("2018-04-30");# End of Florida WY2018

enp.rf=data.frame()
pb=txtProgressBar(1,max=length(dbkey.list$Dbkey),style=3)
for(i in 1:length(dbkey.list$Dbkey)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,dbkey.list$Dbkey[i])
  enp.rf=rbind(tmp,enp.rf)
  setTxtProgressBar(pb,i)
}
nrow(ddply(enp.rf,c("DBKEY"),summarise,N.val=N(DBKEY)))
nrow(ddply(dbkey.list,c("Dbkey"),summarise,N.val=N(Dbkey)))
i
dbkey.list[i,]

rf.dat=ddply(enp.rf,c("Station","Date"),summarise,RF.m=mean(in.to.cm(Data.Value)/100,na.rm=T))
range(rf.dat$RF.m,na.rm=T)
rf.dat$RF.m=with(rf.dat,ifelse(RF.m<0,0,RF.m))
rf.dat$WY=WY(rf.dat$Date)

rf.dat2=ddply(rf.dat,c("Date","WY"),summarise,TRF.m=mean(RF.m,na.rm=T),N.site=N(Station))
rf.dat2$hydro.day=hydro.day(rf.dat2$Date,"FL")
rf.dat2$cum.rf=with(rf.dat2,ave(TRF.m,WY,FUN=function(x)cumsum(x)))
rf.dat2$TRF.cm=rf.dat2$TRF.m*100
rf.dat2$Date.EST=date.fun(rf.dat2$Date)
rf.dat2$DateTime.EST=with(rf.dat2,date.fun(paste(Date.EST,"01:00:00"),form="%Y-%m-%d %H:%M:%S"))#fake Date Time

plot(cum.rf~Date.EST,subset(rf.dat2,WY==2012))
abline(v=date.fun("2011-10-21"))

xlim.val=c(date.fun(wilma.landfall.FL-ddays(10)),date.fun(wilma.landfall.FL+ddays(5)))
plot(TRF.cm~Date,subset(rf.dat2,WY==WY(wilma.landfall.FL)),xlim=xlim.val)
abline(v=wilma.landfall.FL)
subset(rf.dat2,Date==wilma.landfall.FL)
(subset(rf.dat2,Date==wilma.landfall.FL)$TRF.m*gArea(enp))*1e-9

xlim.val=c(date.fun(irma.landfall.FL-ddays(10)),date.fun(irma.landfall.FL+ddays(5)))
plot(TRF.cm~Date,subset(rf.dat2,WY==WY(irma.landfall.FL)),xlim=xlim.val)
abline(v=irma.landfall.FL)
subset(rf.dat2,Date==irma.landfall.FL)
(subset(rf.dat2,Date==irma.landfall.FL)$TRF.m*gArea(enp))*1e-9

quant.sum=ddply(rf.dat2,"hydro.day",summarise,L.Q=quantile(cum.rf,probs=0.10,na.rm=T),U.Q=quantile(cum.rf,probs=0.90,na.rm=T),med=median(cum.rf,na.rm=T))
ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,366);by.x=90;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/3)
#tiff(filename=paste0(plot.path,"FigX_Cumrainfall_bw.tiff"),width=6,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"FigX_Cumrainfall_bw.png"),width=6,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,3.2,0.5,0.5),oma=c(1,1.5,0.5,0.25),mgp=c(3,1,0));
layout(matrix(1:2,2,1,byrow=F),heights=c(0.8,0.2))

plot(TRF.m~hydro.day,rf.dat2,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(quant.sum,shaded.range(hydro.day,L.Q,U.Q,"grey",lty=1))
with(quant.sum,lines(hydro.day,med,lwd=0.5,col="grey"))
with(subset(rf.dat2,WY==WY(wilma.landfall.FL)),lines(hydro.day,cum.rf,lwd=1.75,col="black",lty=2))
with(subset(rf.dat2,WY==WY(wilma.landfall.FL)&hydro.day==hydro.day(wilma.landfall.FL,"FL")),points(hydro.day,cum.rf,pch=19,cex=1.25))
with(subset(rf.dat2,WY==WY(irma.landfall.FL)),lines(hydro.day,cum.rf,col="black",lty=1,lwd=1.75))
with(subset(rf.dat2,WY==WY(irma.landfall.FL)&hydro.day==hydro.day(irma.landfall.FL,"FL")),points(hydro.day,cum.rf,pch=15,cex=1.25))
axis_fun(1,xmaj,xmin,xmaj,1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=2,"Day of water year")
mtext(side=2,line=2.25,"Cumulative Rainfall Total (meters)")

plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
legend.text=c("Interquantile Range (10\u1d57\u02b0 - 90\u1d57\u02b0)",paste("WY",WY(wilma.landfall.FL)),paste("WY",WY(irma.landfall.FL)),NA,"Wilma landfall","Irma landfall")
pt.cols=c(adjustcolor("grey",0.1),NA,NA,NA,"black","black")
legend(0.5,0.5,legend=legend.text,lty=c(0,2,1,0,0,0),col=c("grey","black","black",NA,"black","black"),lwd=c(1,2,2,1,1),pch =c(22,NA,NA,NA,19,15),pt.bg=pt.cols,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()
