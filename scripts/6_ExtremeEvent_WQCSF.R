# BGChem CSF --------------------------------------------------------------
#   Evaluting the affect of hurricane occurance and intensity on water quality 
# across FCE (i.e. SRS and TS)

# Hurricane Track Data
hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,Year.start=min(Year,na.rm=T),WY.start=min(WY,na.rm=T),max.wind.ms=max(Wind*0.44704,na.rm=T),min.pres=min(Pressure,na.rm=T))
hur.dat.sp.pt=SpatialPointsDataFrame(coords=hur.dat[,c("Lon","Lat")],data=hur.dat,proj4string = CRS("+init=epsg:4269"))

#Convert point to line data
hur.dat2=hur.dat
hur.id=ddply(hur.dat,c("Key","Name"),summarise,N.val=length(which(Key!="NA")))
hur.dat2$dLat=with(hur.dat2,ave(Lat,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dLon=with(hur.dat2,ave(Lon,Key,FUN=function(x)c(0,diff(x))))
hur.dat2$dist.m=with(hur.dat2,sqrt((dLat^2)+(dLon^2)));#calculates distance between points

hur.id=ddply(hur.dat2,c("Key","Name"),summarise,N.val=N(Key),max.dist=max(dist.m,na.rm=T))
hur.dat2=subset(hur.dat2,Key%in%subset(hur.id,N.val>1&max.dist>0)$Key);#need greater than one point and some distance between points to create a line
coordinates(hur.dat2)=c("Lon","Lat")
hur.dat2=hur.dat2[order(hur.dat2$Key,hur.dat2$DateTime),]
hur.dat2$WY=WY(hur.dat2$DateTime)
path=sp::split(hur.dat2,hur.dat2$Key)

##Convert spatial points to spatial lines for each hurricane
sp_lines=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
pb=txtProgressBar(1,max=length(path),style=3)
for(i in 2:length(path)){
  tmp=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
  #tmp=SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269"))
  sp_lines=rbind(sp_lines,tmp)
  setTxtProgressBar(pb,i)
}
chk=data.frame(gIsValid(sp_lines,byid=T));#checks
chk$Key=rownames(chk)
colnames(chk)=c("geo.valid","Key")
subset(chk,geo.valid=="FALSE")

hur.tck=spTransform(sp_lines,utm17)
hur.tck=merge(hur.tck,hur.year,by.x=c("Key","Name"),by.y=c("Key","Name"))

#ENP boundary and buffer
enp.buffer=gBuffer(enp,width=200*1000);# 200 km buffer around ENP
enp.buffer=SpatialPolygonsDataFrame(enp.buffer,data.frame(row.names = "buffer",width.km=200,area="ENP"))

enp.hurr=over(enp.buffer,hur.tck,returnList = T,byid=T)[[1]];#select only hurricanes that cross the 200 km buffer. 
enp.hurr=enp.hurr[order(enp.hurr$WY.start,enp.hurr$Key),]

#hur.dat.enp=subset(hur.dat,Key%in%enp.hurr$Key&WY%in%seq(2005,2018,1))
hur.dat.enp=subset(hur.dat,Key%in%enp.hurr$Key)

## Accumulated Cyclone Energy (ACE)
max.wind=ddply(hur.dat.enp,c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
ENP.ACE=ddply(max.wind,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))
ENP.ACE=merge(ENP.ACE,data.frame(WY=seq(1853,2018,1),fill=1),"WY",all.y=T)
ENP.ACE$ACE=with(ENP.ACE,ifelse(is.na(ACE)==T,0,ACE))

ENP.ACE=rbind(ENP.ACE,data.frame(WY=c(2008,2010,2012,2015),ACE=0))
ENP.ACE=ENP.ACE[order(ENP.ACE$WY),]

max.wind.all=ddply(hur.dat,c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
ACE.all=ddply(max.wind.all,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))

plot(ACE~WY,ACE.all,ylim=c(0,350))
with(ENP.ACE,points(WY,ACE,pch=21,bg="red"))
abline(v=WY(c(wilma.landfall.FL,irma.landfall.FL)))

xlim.val=c(2005,2018);by.x=4;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(0,350);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
#tiff(filename=paste0(plot.path,"FigX_HurricaneACE.tiff"),width=6,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,2.5,0,0.25),oma=c(1.5,1.5,1.5,0.5),mgp=c(3,1,0));
layout(matrix(c(1:2),2,1,byrow=F),heights=c(0.8,0.2))

plot(ACE~WY,ACE.all,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
polygon(c(2004,2004,2019,2019),c(66,111,111,66),col=adjustcolor("grey",0.25),border="grey")#density=20)
with(ACE.all,pt_line(WY,ACE,2,"black",1,21,"grey",cex=1.5))
with(ENP.ACE,pt_line(WY,ACE,2,"indianred1",1,21,"indianred1",cex=1.5))
axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,line=2,"Water Year")
mtext(side=2,line=2.5,"Accumulated Cyclone Energy (x10\u2074 knots\u00B2 Yr\u207B\u00B9)")

plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
legend.text=c("Atlantic Basin ACE","ACE w/in 200-km radius of ENP","Atlantic Basin\nNear-normal season (66 > ACE < 111)")
legend(0.5,0.25,legend=legend.text,lty=c(2,2,NA),col=c("black","indianred1","grey"),lwd=c(1,1,NA),pch =c(21,21,22),pt.bg=c("grey","indianred1",adjustcolor("grey",0.25)),pt.cex=1.25,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col = "white")
legend(0.5,0.25,legend=legend.text,lty=c(NA,NA,NA),col=c("black","black"),lwd=c(1,1,NA),pch =c(21,21,NA),pt.bg=c("grey","indianred1","white"),pt.cex=1.25,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

# Water Quality data ------------------------------------------------------

## FCE data (grab data)
fce.grahl=read.csv(paste0(data.path,"/WQ/LT_ND_Grahl_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.losada=read.csv(paste0(data.path,"/WQ/LT_ND_Losada_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rondeau=read.csv(paste0(data.path,"/WQ/LT_ND_Rondeau_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rubio=read.csv(paste0(data.path,"/WQ/LT_ND_Rubio_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))

names(fce.grahl)
names(fce.losada)
names(fce.rondeau)
names(fce.rubio)
fce.grahl=rename(fce.grahl,c("NandN"="N.N"));#renaming NOx column to be consistent with other datasets
fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio)
fce.wq$Date.EST=date.fun(as.character(fce.wq$DATE),form="%Y-%m-%d")
#subset(fce.wq,is.na(Date.EST)==T); #sanity check

unique(fce.wq$SITENAME)
fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Region=c(rep("SRS",8),rep("TS",8)))

fce.wq=merge(fce.wq,fce.wq.sites,"SITENAME")
unique(fce.wq$SITENAME)

fce.wq$TP=round((fce.wq$TP*P.mw),4);#convert uM concentration to ug/L
fce.wq$TN=round((fce.wq$TN*N.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$TOC=round((fce.wq$TOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$DOC=round((fce.wq$DOC*C.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq$SRP=round((fce.wq$SRP*P.mw),4);#convert uM concentration to ug/L
fce.wq$SRP=with(fce.wq,ifelse(SRP<=0.3,0.3097,SRP))
fce.wq$N.N=with(fce.wq,ifelse(N.N<=0,0.01,N.N))
fce.wq$DIN=with(fce.wq,round(((N.N+NH4)*N.mw)*0.001,4));#convert uM concentration to mg/L

vars3=c("Station.ID","Region","Date.EST","TP","TN","TOC","DOC","SRP","DIN")
fce.wq=fce.wq[,vars3]
head(fce.wq)

grab.wq=fce.wq

# spatial data
fce.sp.site.list=c(paste0("SRS-",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/Ph-",c("1a","1b","2a","2b",2,3,"6b","6a","7a","7b")))
fce.sites2=subset(fce.sites,SITE%in%fce.sp.site.list)
fce.sites2=merge(fce.sites2,data.frame(SITE=fce.sp.site.list,STATION=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b","2a","2b",2,3,"6b","6a","7a","7b"))),Region=c(rep("SRS",8),rep("TS",10))),"SITE")
fce.sites2$datasource="FCE"
wq.sites.sp=fce.sites2[,c("STATION","Region","datasource")]

TS2.locdat=data.frame(STATION="TS/PH2",Region="TS",datasource="FCE",LONG=-80.607,LAT=25.404)
TS2=spTransform(SpatialPointsDataFrame(coords=TS2.locdat[,c("LONG","LAT")],data=TS2.locdat[,c("STATION","Region","datasource")],proj4string=CRS("+init=epsg:4326")),utm17)
wq.sites.sp=union(wq.sites.sp,TS2)

tmap_mode("view")
tm_shape(fce.sites2)+tm_dots()+
  tm_shape(TS2)+tm_dots(col="red")
tm_shape(wq.sites.sp)+tm_dots(col="datasource",size=0.1,palette=c("blue","green"))

##
## Euclidean Distance from GOM/FB
srs.struct=subset(sfwmd.mon.q,SITE%in%c(paste0("S12",LETTERS[1:4]),"S333"))
SRSr=raster(xmn=extent(enp)[1],xmx=extent(enp)[2],ymn=extent(enp)[3],ymx=extent(enp)[4],crs=utm17)
SRSr=setValues(SRSr,0)
SRSr=mask(SRSr,srs.struct)
SRSrD=distance(SRSr)
plot(SRSrD);plot(srs.struct,add=T)
wq.sites.sp.srs=subset(wq.sites.sp,Region=="SRS")
wq.sites.sp.srs$EuDist.m=raster::extract(SRSrD,wq.sites.sp.srs)
wq.sites.sp.srs$EuDist.frac=wq.sites.sp.srs$EuDist.m/max(wq.sites.sp.srs$EuDist.m)

ts.struct=subset(sfwmd.mon.q,SITE%in%c("S332D","S18C","S199"))
TSr=raster(xmn=extent(enp)[1],xmx=extent(enp)[2],ymn=extent(enp)[3],ymx=extent(enp)[4],crs=utm17)
TSr=setValues(TSr,0)
TSr=mask(TSr,ts.struct)
TSrD=distance(TSr)
plot(TSrD);plot(ts.struct,add=T)
wq.sites.sp.ts=subset(wq.sites.sp,Region=="TS")
wq.sites.sp.ts$EuDist.m=raster::extract(TSrD,wq.sites.sp.ts)
wq.sites.sp.ts$EuDist.frac=wq.sites.sp.ts$EuDist.m/max(wq.sites.sp.ts$EuDist.m)

wq.sites.sp=rbind(wq.sites.sp.srs,wq.sites.sp.ts)
#head(wq.sites.sp@data)

# Merging WQ and Euclidean distance data
grab.wq2=merge(grab.wq,wq.sites.sp@data[,c("STATION","datasource","EuDist.frac")],by.x="Station.ID",by.y="STATION",all.x=T)
head(grab.wq2)
grab.wq2$month=format(grab.wq2$Date.EST,"%m")
grab.wq2$CY=as.numeric(format(grab.wq2$Date.EST,"%Y"))
grab.wq2$WY=WY(grab.wq2$Date.EST)
grab.wq2$dec.year=decimal_date(grab.wq2$Date.EST)
grab.wq2=grab.wq2[order(grab.wq2$Station.ID,grab.wq2$Date.EST),]
grab.wq2$SRP=with(grab.wq2,ifelse(SRP<=0,"2",SRP))
#grab.wq2=merge(grab.wq2,fce.wq.sites[,c("Station.ID","Region")],"Station.ID")

## Annual trend analysis
site.val=c(paste0("SRS",c("1d",2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))

WQ.WY=ddply(grab.wq2,c("Station.ID","Region","datasource","EuDist.frac","WY"),summarise,mean.TP=mean(TP,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),N.TN=N(TN),mean.DOC=mean(DOC,na.rm=T),N.DOC=N(DOC),mean.SRP=mean(SRP,na.rm=T),N.SRP=N(SRP),mean.DIN=mean(DIN,na.rm=T),N.DIN=N(DIN))
WQ.WY$TP.scn=with(WQ.WY,ifelse(N.TP<=4,NA,mean.TP))
WQ.WY$TN.scn=with(WQ.WY,ifelse(N.TN<=4,NA,mean.TN))
WQ.WY$DOC.scn=with(WQ.WY,ifelse(N.DOC<=4,NA,mean.DOC))
WQ.WY$SRP.scn=with(WQ.WY,ifelse(N.SRP<=4,NA,mean.SRP))
WQ.WY$DIN.scn=with(WQ.WY,ifelse(N.DIN<=4,NA,mean.DIN))
WQ.WY$dec.yr=with(WQ.WY,WY+0.328);#decimal date for Apirl 30th 
WQ.WY$WY.f=as.factor(WQ.WY$WY)
WQ.WY$WY.f2=paste0("WY",WQ.WY$WY)
WQ.WY$Hurr=as.factor(with(WQ.WY,ifelse(WY%in%c(2006,2018),1,0)))
WQ.WY=subset(WQ.WY,Station.ID%in%site.val)

plot(TP.scn~WY,subset(WQ.WY,Station.ID==site.val[1]))
plot(TP.scn~WY,subset(WQ.WY,Station.ID==site.val[2]))
plot(TP.scn~WY,subset(WQ.WY,Station.ID==site.val[3]))
plot(TP.scn~WY,subset(WQ.WY,Station.ID==site.val[4]))

## Time Series plot
unique(WQ.WY$Station.ID)
site.val=c(paste0("SRS",c("1d",2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))
WQ.WY$Station.ID=factor(WQ.WY$Station.ID,level=site.val)
#tiff(filename=paste0(plot.path,"FigX_WQPlots_v2.tiff"),width=9,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))
layout(matrix(1:48,6,8,byrow=F),widths=c(1,1,0.3,1,1,0.30,1,1))

xlim.val=c(2001.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(1,110);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TP.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),TP,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,mean.TP,pch=21,bg=ifelse(is.na(TP.scn)==T,NA,"indianred1"),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1993.5,y=ylim.val[1],"Total Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0.1,5);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TN.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),TN,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,mean.TN,pch=22,bg=ifelse(is.na(TN.scn)==T,NA,"dodgerblue1"),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1993.5,y=ylim.val[1],"Total Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(1,55);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(DOC.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),DOC,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,mean.DOC,pch=23,bg=ifelse(is.na(DOC.scn)==T,NA,"darkgoldenrod1"),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1993.5,y=ylim.val[2],"Dissolved Organic Carbon (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,line=0.5,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean TP", "Annual Mean TN","Annual Mean DOC", "Grab Samples")
pt.cols=c("indianred1","dodgerblue1","darkgoldenrod1",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,23,19),col=c("black","black","black",pt.cols[4]),lwd=0.2,lty=NA,pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

unique(WQ.WY$Station.ID)
site.val=c(paste0("SRS",c("1d",2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))
WQ.WY$Station.ID=factor(WQ.WY$Station.ID,level=site.val)
#tiff(filename=paste0(plot.path,"FigX_WQPlots_inorganic.tiff"),width=7,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,1))
#layout(matrix(1:12,6,2,byrow=F))
layout(matrix(1:30,6,5,byrow=F),widths=c(1,1,0.3,1,1))

xlim.val=c(2001.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(0.1,65);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TP.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),SRP,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,mean.SRP,pch=21,bg=ifelse(is.na(SRP.scn)==T,NA,"sienna1"),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995.5,y=ylim.val[2],"Soluable Reactive Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0.01,3);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TN.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),DIN,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,mean.DIN,pch=22,bg=ifelse(is.na(DIN.scn)==T,NA,"darkcyan"),cex=1.25,lwd=0.2))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==5){text(x=1995,y=ylim.val[2],"Dissolved Inorganic Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean SRP", "Annual Mean DIN", "Grab Samples")
pt.cols=c("sienna1","darkcyan",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,19),col=c("black","black",pt.cols[3]),lwd=0.2,lty=NA,pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

#tiff(filename=paste0(plot.path,"FigX_WQPlots.tiff"),width=6.5,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2.75,0.5,1))
#layout(matrix(1:12,6,2,byrow=F))
layout(matrix(1:30,6,5,byrow=F),widths=c(1,1,0.3,1,1))

xlim.val=c(2004.8,2018.5);by.x=5;xmaj=seq(round(xlim.val[1]),xlim.val[2],by.x);xmin=seq(round(xlim.val[1]),xlim.val[2],by.x/by.x)
ylim.val=c(1,110);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TP.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),TP,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,TP.scn,pch=21,bg="indianred1",cex=1.5))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1999,y=ylim.val[1],"Total Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
for(j in 1:7){plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)}

ylim.val=c(0.1,5);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:length(site.val)){
  plot(TN.scn~WY,WQ.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
  with(subset(grab.wq2,Station.ID==site.val[i]),points(jitter(WY,1),TN,pch=19,col=adjustcolor("grey",0.75)))
  with(subset(WQ.WY,Station.ID==site.val[i]),points(WY,TN.scn,pch=22,bg="dodgerblue1",cex=1.5))
  if(i==6|i==11){axis_fun(1,line=-0.5,xmaj,xmin,xmaj)}else(axis_fun(1,xmaj,xmin,NA))
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,site.val[i],cex=0.8)
  if(i==4){text(x=1999.75,y=ylim.val[1],"Dissolved Inorganic Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean TP", "Annual Mean TN", "Grab Samples")
pt.cols=c("indianred1","dodgerblue1",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,19),col=c("black","black",pt.cols[3]),pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()


#Climate CSF
WQ.WY2=subset(WQ.WY,WY%in%seq(2002,2018,1))
WQ.WY2=merge(WQ.WY2,ENP.ACE,'WY')

WQ.WY2.TP=subset(WQ.WY2,is.na(TP.scn)==F)
subset(WQ.WY2.TP,is.na(Station.ID)==T)

WQ.WY2.TN=subset(WQ.WY2,is.na(TN.scn)==F)
subset(WQ.WY2.TN,is.na(Station.ID)==T)

WQ.WY2.DOC=subset(WQ.WY2,is.na(DOC.scn)==F)
subset(WQ.WY2.DOC,is.na(Station.ID)==T)

#######MODEL SELECTION EXERCISE: NONLINEARITY ----------
#nonlinear Climate Sensitivity Functions (CSFs)
#Does TP vary nonlinearly with the hurricane intensity index (ACE)?

#Nonlinear responses to climate drivers are a way to assess
#whether TP will be sensitive (or not) to increasing VARIABILITY in Hurricane intensity
#Nonlinearities indicate sensitivity to change in variance of climate predictor
#Linear suggests no sensitivity to to change in variance of climate predictor

#null, 1-linear, 2-quadratic, 3-cubic; AR = autoregressive 1, .y = adds random effect of year

## FCE wide models
# TP
# Non-linear mixed model
mnull.TP<-lme(log(TP.scn)~1,data=WQ.WY2.TP,random=~1|Station.ID,method="ML")
m1.TP=lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=~1|Station.ID,method="ML")
m2.TP<-update(m1.TP, .~. +I(ACE^2))
m3.TP<-update(m2.TP, .~. +I(ACE^3))

AICc.TP=AICc(mnull.TP,m1.TP, m2.TP, m3.TP)
subset(AICc.TP,AICc==min(AICc.TP$AICc))
model.sel(mnull.TP,m1.TP, m2.TP, m3.TP)

# variance-covariance matrix first-order autoregressive var-cov matrix
mnullAR.TP<-lme(log(TP.scn)~1,data=WQ.WY2.TP,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m1AR.TP<-lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TP<-update(m1AR.TP, .~. +I(ACE^2))
m3AR.TP<-update(m2AR.TP, .~. +I(ACE^3))

AICc.TP=AICc(mnullAR.TP,m1AR.TP, m2AR.TP, m3AR.TP)
subset(AICc.TP,AICc==min(AICc.TP$AICc))
model.sel(mnullAR.TP,m1AR.TP,m2AR.TP,m3AR.TP)

AICc(mnull.TP,m1.TP, m2.TP, m3.TP,mnullAR.TP,m1AR.TP, m2AR.TP, m3AR.TP)
model.sel(mnull.TP,m1.TP, m2.TP, m3.TP,mnullAR.TP,m1AR.TP,m2AR.TP,m3AR.TP)

anova(m2AR.TP,m3AR.TP);#no significant difference between models, use the simplest. 

mod.TP=visreg(m2AR.TP,"ACE",trans=exp,type="conditional")

Anova(m2AR.TP,type=3)
summary(m2AR.TP)
rsquared(m2AR.TP)

#double check model assumptions
plot(m2AR.TP); #residual versus fitted (qualitative)
qqnorm(residuals(m2AR.TP));qqline(residuals(m2AR.TP)); #normality of residuals (qualitative)
shapiro.test(residuals(m2AR.TP)); #normality of residuals


#m2AR2.TP=update(m2AR.TP,correlation=corARMA(form=~WY,p=2))
#visreg(m2AR2.TP,"ACE",trans=exp,type="conditional")
#rsquared(m2AR2.TP)
#AICc(m2AR2.TP)
#AR2 is not improved

# CSF models with .y = year as a random effect 
# added this previously but decided to not include in full analysis
#mnull.WY.TP<-lme(log(TP.scn)~1,data=WQ.WY2.TP,random=list(~1|Station.ID,~1|WY.f),method="ML")
#m1.WY.TP=lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=list(~1|Station.ID,~1|WY.f),method="ML")
#m2.WY.TP<-update(m1.WY.TP, .~. +I(ACE^2))
#m3.WY.TP<-update(m2.WY.TP, .~. +I(ACE^3))
#m1AR.WY.TP<-lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=list(~1|Station.ID,~1|WY.f),correlation=corAR1(form=~WY),method="ML")
#m2AR.WY.TP<-update(m1AR.WY.TP, .~. +I(ACE^2))
#m3AR.WY.TP<-update(m2AR.WY.TP, .~. +I(ACE^3))

# TN
# Non-linear mixed model
mnull.TN<-lme(log(TN.scn)~1,data=WQ.WY2.TN,random=~1|Station.ID,method="ML")
m1.TN=lme(log(TN.scn)~ACE,data=WQ.WY2.TN,random=~1|Station.ID,method="ML")
m2.TN<-update(m1.TN, .~. +I(ACE^2))
m3.TN<-update(m2.TN, .~. +I(ACE^3))

AICc.TN=AICc(mnull.TN,m1.TN, m2.TN, m3.TN)
subset(AICc.TN,AICc==min(AICc.TN$AICc))
model.sel(mnull.TN,m1.TN, m2.TN, m3.TN)

# variance-covariance matrix first-order autoregressive var-cov matrix
mnullAR.TN<-lme(log(TN.scn)~1,data=WQ.WY2.TN,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m1AR.TN<-lme(log(TN.scn)~ACE,data=WQ.WY2.TN,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TN<-update(m1AR.TN, .~. +I(ACE^2))
m3AR.TN<-update(m2AR.TN, .~. +I(ACE^3))

AICc.TN=AICc(mnullAR.TN,m1AR.TN, m2AR.TN, m3AR.TN)
subset(AICc.TN,AICc==min(AICc.TN$AICc))
model.sel(mnullAR.TN,m1AR.TN,m2AR.TN,m3AR.TN)

AICc(mnull.TN,m1.TN, m2.TN, m3.TN,mnullAR.TN,m1AR.TN, m2AR.TN, m3AR.TN)
model.sel(mnull.TN,m1.TN, m2.TN, m3.TN,mnullAR.TN,m1AR.TN,m2AR.TN,m3AR.TN)

anova(m1AR.TN,m2AR.TN);#no difference between model, therefore go for the simpliest. 
anova(m1AR.TN,mnullAR.TN)
mod.TN=visreg(m1AR.TN,"ACE",trans=exp,type="conditional")
visreg(m2AR.TN,"ACE",trans=exp,type="conditional")

Anova(m1AR.TN,type=3)
summary(m1AR.TN)
rsquared(m1AR.TN)

#double check model assumptions
plot(m1AR.TN); #residual versus fitted (qualitative)
qqnorm(residuals(m1AR.TN));qqline(residuals(m1AR.TN)); #normality of residuals (qualitative)
shapiro.test(residuals(m1AR.TN)); #normality of residuals

# DOC
# Non-linear mixed model
mnull.DOC<-lme(log(DOC.scn)~1,data=WQ.WY2.DOC,random=~1|Station.ID,method="ML")
m1.DOC=lme(log(DOC.scn)~ACE,data=WQ.WY2.DOC,random=~1|Station.ID,method="ML")
m2.DOC<-update(m1.DOC, .~. +I(ACE^2))
m3.DOC<-update(m2.DOC, .~. +I(ACE^3))

AICc.DOC=AICc(mnull.DOC,m1.DOC, m2.DOC, m3.DOC)
subset(AICc.DOC,AICc==min(AICc.DOC$AICc))
model.sel(mnull.DOC,m1.DOC, m2.DOC, m3.DOC)

# variance-covariance matrix first-order autoregressive var-cov matrix
mnullAR.DOC<-lme(log(DOC.scn)~1,data=WQ.WY2.DOC,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m1AR.DOC<-lme(log(DOC.scn)~ACE,data=WQ.WY2.DOC,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.DOC<-update(m1AR.DOC, .~. +I(ACE^2))
m3AR.DOC<-update(m2AR.DOC, .~. +I(ACE^3))

AICc.DOC=AICc(mnullAR.DOC,m1AR.DOC, m2AR.DOC, m3AR.DOC)
subset(AICc.DOC,AICc==min(AICc.DOC$AICc))
model.sel(mnullAR.DOC,m1AR.DOC,m2AR.DOC,m3AR.DOC)

AICc(mnull.DOC,m1.DOC, m2.DOC, m3.DOC,mnullAR.DOC,m1AR.DOC, m2AR.DOC, m3AR.DOC)
model.sel(mnull.DOC,m1.DOC, m2.DOC, m3.DOC,mnullAR.DOC,m1AR.DOC,m2AR.DOC,m3AR.DOC)
anova(m1.DOC,m1AR.DOC)

mod.DOC=visreg(m1.DOC,"ACE",trans=exp,type="conditional")
visreg(m1AR.DOC,"ACE",trans=exp,type="conditional")

Anova(m1.DOC,type=3)
summary(m1.DOC)
rsquared(m1.DOC)
rsquared(m1AR.DOC)


#tiff(filename=paste0(plot.path,"FigX_CSF_all.tiff"),width=4,height=6,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.25,1.75,1,0.25),mar=c(1.5,2,0.5,1))
layout(matrix(1:4,4,1,byrow=F),heights=c(rep(1,3),0.25))
xlim.val=c(0,180);by.x=25;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
axis.lab.cex=0.8
{
ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)#;ymaj=c(ylim.val[1],log.scale.fun(ylim.val[1],ylim.val[2],"major"));ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")#
plot(mean.TP~ACE,WQ.WY2,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(v=xmaj,h=ymaj,lty=3,col="grey")
with(WQ.WY2,points(ACE,mean.TP,pch=ifelse(Region=="SRS",21,22),bg=ifelse(Region=="SRS","indianred1","dodgerblue1"),lwd=0.2))
with(mod.TP$fit,shaded.range(ACE,visregLwr,visregUpr,"grey",lty=0,col.adj=0.5))
with(mod.TP$fit,lines(ACE,visregFit,lwd=1.25,lty=2))
axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Annual Mean TP (\u03BCg L\u207B\u00B9)",cex=axis.lab.cex)

ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)#;ymaj=c(ylim.val[1],log.scale.fun(ylim.val[1],ylim.val[2],"major"));ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")#
plot(mean.TN~ACE,WQ.WY2,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(v=xmaj,h=ymaj,lty=3,col="grey")
with(WQ.WY2,points(ACE,mean.TN,pch=ifelse(Region=="SRS",21,22),bg=ifelse(Region=="SRS","indianred1","dodgerblue1"),lwd=0.2))
with(mod.TN$fit,shaded.range(ACE,visregLwr,visregUpr,"grey",lty=0,col.adj=0.5))
with(mod.TN$fit,lines(ACE,visregFit,lwd=1.25,lty=2))
axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Annual Mean TN (mg L\u207B\u00B9)",cex=axis.lab.cex)

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)#;ymaj=c(ylim.val[1],log.scale.fun(ylim.val[1],ylim.val[2],"major"));ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")#
plot(mean.DOC~ACE,WQ.WY2,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(v=xmaj,h=ymaj,lty=3,col="grey")
with(WQ.WY2,points(ACE,mean.DOC,pch=ifelse(Region=="SRS",21,22),bg=ifelse(Region=="SRS","indianred1","dodgerblue1"),lwd=0.2))
with(mod.DOC$fit,shaded.range(ACE,visregLwr,visregUpr,"grey",lty=0,col.adj=0.5))
with(mod.DOC$fit,lines(ACE,visregFit,lwd=1.25,lty=2))
axis_fun(1,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Annual Mean DOC (mg L\u207B\u00B9)",cex=axis.lab.cex)
mtext(side=1,line=2,"Accumulated Cyclone Energy (x10\u2074 knots\u00B2 Yr\u207B\u00B9)",cex=axis.lab.cex)
}
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Shark River Slough","Taylor Slough","Non-linear Mixed Model \u00B1 95% CI")
pt.cols=c("indianred1","dodgerblue1",adjustcolor("grey",0.5))
legend(0.5,0.25,legend=legend.text,pch =c(21,22,22),col=c("black","black","grey"),lwd=0.2,lty=NA,
       pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()


##
##
##
##
##
##
##
##
##
##

# TP
# Shark River Slough
# Non-linear mixed model
mnull.TP.SRS<-lme(log(TP.scn)~1,data=subset(WQ.WY2.TP,Region=="SRS"),random=~1|Station.ID,method="ML")
m1.TP.SRS=lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="SRS"),random=~1|Station.ID,method="ML")
m2.TP.SRS<-update(m1.TP.SRS, .~. +I(ACE^2))
m3.TP.SRS<-update(m2.TP.SRS, .~. +I(ACE^3))

AICc.TP.SRS=AICc(mnull.TP.SRS,m1.TP.SRS, m2.TP.SRS, m3.TP.SRS)
subset(AICc.TP.SRS,AICc==min(AICc.TP.SRS$AICc))
model.sel(mnull.TP.SRS,m1.TP.SRS, m2.TP.SRS, m3.TP.SRS)

# variance-covariance matrix first-order autoregressive var-cov matrix
mnullAR.TP.SRS<-lme(log(TP.scn)~1,data=subset(WQ.WY2.TP,Region=="SRS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m1AR.TP.SRS<-lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="SRS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TP.SRS<-update(m1AR.TP.SRS, .~. +I(ACE^2))
m3AR.TP.SRS<-update(m2AR.TP.SRS, .~. +I(ACE^3))

AICc.TP.SRS=AICc(mnullAR.TP.SRS,m1AR.TP.SRS, m2AR.TP.SRS, m3AR.TP.SRS)
subset(AICc.TP.SRS,AICc==min(AICc.TP.SRS$AICc))
model.sel(mnullAR.TP.SRS,m1AR.TP.SRS,m2AR.TP.SRS,m3AR.TP.SRS)

model.sel(mnull.TP.SRS,m1.TP.SRS, m2.TP.SRS, m3.TP.SRS,m1AR.TP.SRS,m2AR.TP.SRS,m3AR.TP.SRS)

anova(m2AR.TP.SRS,m3AR.TP.SRS);#no significant difference between models, use the simplest. 

mod.TP=visreg(m2AR.TP.SRS,"ACE",trans=exp,type="conditional")

Anova(m2AR.TP.SRS,type=3)
summary(m2AR.TP.SRS)
rsquared(m2AR.TP.SRS)

#double check model assumptions
plot(m2AR.TP.SRS); #residual versus fitted (qualitative)
qqnorm(residuals(m2AR.TP.SRS));qqline(residuals(m2AR.TP.SRS)); #normality of residuals (qualitative)
shapiro.test(residuals(m2AR.TP.SRS)); #normality of residuals
## Does not fit assumptions of the test 

# Taylor Slough
# Non-linear mixed model
mnull.TP.TS<-lme(log(TP.scn)~1,data=subset(WQ.WY2.TP,Region=="TS"),random=~1|Station.ID,method="ML")
m1.TP.TS=lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="TS"),random=~1|Station.ID,method="ML")
m2.TP.TS<-update(m1.TP.TS, .~. +I(ACE^2))
m3.TP.TS<-update(m2.TP.TS, .~. +I(ACE^3))

AICc.TP.TS=AICc(mnull.TP.TS,m1.TP.TS, m2.TP.TS, m3.TP.TS)
subset(AICc.TP.TS,AICc==min(AICc.TP.TS$AICc))
model.sel(mnull.TP.TS,m1.TP.TS, m2.TP.TS, m3.TP.TS)

# variance-covariance matrix first-order autoregressive var-cov matrix
mnullAR.TP.TS<-lme(log(TP.scn)~1,data=subset(WQ.WY2.TP,Region=="TS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m1AR.TP.TS<-lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="TS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TP.TS<-update(m1AR.TP.TS, .~. +I(ACE^2))
m3AR.TP.TS<-update(m2AR.TP.TS, .~. +I(ACE^3))

AICc.TP.TS=AICc(mnullAR.TP.TS,m1AR.TP.TS, m2AR.TP.TS, m3AR.TP.TS)
subset(AICc.TP.TS,AICc==min(AICc.TP.TS$AICc))
model.sel(mnullAR.TP.TS,m1AR.TP.TS,m2AR.TP.TS,m3AR.TP.TS)

model.sel(mnull.TP.TS,m1.TP.TS, m2.TP.TS, m3.TP.TS,m1AR.TP.TS,m2AR.TP.TS,m3AR.TP.TS)

anova(m2AR.TP.TS,m2.TP.TS);#no significant difference between models, use the simplest. 

mod.TP=visreg(m2.TP.TS,"ACE",trans=exp,type="conditional")

Anova(m2.TP.TS,type=3)
summary(m2.TP.TS)
rsquared(m2.TP.TS)

#double check model assumptions
plot(m2.TP.TS); #residual versus fitted (qualitative)
qqnorm(residuals(m2.TP.TS));qqline(residuals(m2.TP.TS)); #normality of residuals (qualitative)
shapiro.test(residuals(m2.TP.TS)); #normality of residuals
#Does fit assumptions of the test


##
##
##
##
##
##
##
##
##
##
