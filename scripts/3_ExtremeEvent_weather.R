# 

# Weather Data ------------------------------------------------------------
wx.dat2005=read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=vcaf1h2005.txt.gz&dir=data/historical/stdmet/",sep="",header=T,na.strings=c("99","999","9999.0"))
wx.dat2017=read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=vcaf1h2017.txt.gz&dir=data/historical/stdmet/",sep="",header=F,na.strings=c("99","999","9999.0"))
#wx.dat2005=read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=lonf1h2005.txt.gz&dir=data/historical/stdmet/",sep="",header=T,na.strings=c("99","999","9999.0"))
#wx.dat2017=read.table("https://www.ndbc.noaa.gov/view_text_file.php?filename=lonf1h2017.txt.gz&dir=data/historical/stdmet/",sep="",header=F,na.strings=c("99","999","9999.0"))
colnames(wx.dat2017)=names(wx.dat2005)
wx.dat=rbind(wx.dat2005,wx.dat2017)
wx.dat$DateTime.EST=with(wx.dat,date.fun(paste(YYYY,"-",MM,"-",DD," ",hh,":",ifelse(mm==0,"00",mm),":","00",sep=""),form="%Y-%m-%d %H:%M:%S"))
wx.dat$Date.EST=date.fun(wx.dat$DateTime.EST)
wx.dat$WSPD[wx.dat$WSPD>=99]=NA
wx.dat$ATMP[wx.dat$ATMP>=999.0]=NA
wx.dat$WTMP[wx.dat$WTMP>=999.0]=NA
wx.dat$BP.mmHg=hPa.to.mmHg(wx.dat$BAR)

plot(WSPD~DateTime.EST,subset(wx.dat,YYYY==2017))



landfall=c(wilma.landfall.FL,irma.landfall.FL)
storm.name=c("Wilma","Irma")
#tiff(filename=paste0(plot.path,"Fig2_Weatherinfo.tiff"),width=6.25,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"Fig2_Weatherinfo.png"),width=6.25,height=6,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2.5,0,2),oma=c(2.5,1.5,1.5,2),mgp=c(3,1,0));
layout(matrix(c(1:6),3,2,byrow=F))
y1.lab=2.5
for(i in 1:2){
  hur.poly.time=c(rep(date.fun(paste(landfall[i],"01:00:00"),form="%Y-%m-%d %H:%M:%S"),2),rep(date.fun(paste(landfall[i],"24:00:00"),form="%Y-%m-%d %H:%M:%S"),2))
  period.wx=seq(date.fun(landfall[i]-ddays(10)),date.fun(landfall[i]+ddays(10)),"1 days")
  xlim.val=range(period.wx);xmaj=seq(xlim.val[1],xlim.val[2],"5 days");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
  
  ylim.val=list(c(0,10),c(0,35))[[i]];by.y=list(2,10)[[i]];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(TRF.cm~DateTime.EST,rf.dat2,type="n",xlim=xlim.val,ylim=ylim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,yaxs="i")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  #polygon(x=hur.poly.time,y=c(-1,15,15,-1),border=F,col=adjustcolor("grey90",0.25))
  for(j in 1:nrow(rf.dat2)){
    xx=with(rf.dat2[j,],c(rep(date.fun(paste(Date.EST,"01:00:00"),form="%Y-%m-%d %H:%M:%S"),2),rep(date.fun(paste(Date.EST,"23:00:00"),form="%Y-%m-%d %H:%M:%S"),2)))
    yy=with(rf.dat2[j,],c(0,TRF.cm,TRF.cm,0))
    polygon(xx,yy,col="grey")
  }
  axis_fun(1,xmaj,xmin,NA,1);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
  if(i==1){mtext(side=2,line=y1.lab,"Daily Total Rainfall (cm)")}
  mtext(side=3,paste(storm.name[i]," (",format(landfall[i],"%Y"),")",sep=""))
  
  ylim.val=c(720,770);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(BP.mmHg~DateTime.EST,wx.dat,xlim=xlim.val,ylim=ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  polygon(x=hur.poly.time,y=c(500,800,800,500),border=F,col=adjustcolor("grey90",0.5))
  with(wx.dat,lines(DateTime.EST,BP.mmHg,lwd=1.25))
  axis_fun(2,ymaj,ymin,ymaj,1)
  if(i==1){mtext(side=2,line=y1.lab,"Bar. Press. (mmHg)")}
  ylim.val=list(c(10,35),c(20,40))[[i]];by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  par(new=T);plot(ATMP~DateTime.EST,wx.dat,xlim=xlim.val,ylim=ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
  with(wx.dat,lines(DateTime.EST,ATMP,lty=2,lwd=1.25,col="grey"))
  with(wx.dat,lines(DateTime.EST,WTMP,lty=1,lwd=1.25,col="dodgerblue1"))
  axis_fun(1,xmaj,xmin,NA,1)
  axis_fun(4,ymaj,ymin,ymaj,1)
  if(i==2){mtext(side=4,line=2,"Temperature (\u00B0C)")}
  if(i==1){legend("bottomleft",legend=c("Barometric Pressure","Air Temp","Water Temp"),lty=c(1,2,1),col=c("black","grey","dodgerblue1"),lwd=1.25,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)}
  
  ylim.val=c(-20,25);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(WSPD~DateTime.EST,wx.dat,xlim=xlim.val,ylim=ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  polygon(x=hur.poly.time,y=c(-50,50,50,-50),border=F,col=adjustcolor("grey90",0.5))
  stickplot.dat.arrows(DateTime.EST,WSPD,WD,wx.dat,col="grey50")
  axis_fun(2,ymaj,ymin,abs(ymaj),1)
  if(i==1){mtext(side=2,line=y1.lab,"Wind Velocity (m s\u207B\u00B9)")}
  axis_fun(1,xmaj,xmin,format(xmaj,"%m-%d"),1)
}
mtext(side=1,line=1,"Date (Month-Day)",outer=T)
dev.off()
