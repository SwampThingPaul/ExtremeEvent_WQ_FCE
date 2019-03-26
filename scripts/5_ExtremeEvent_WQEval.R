# Disturbance CDF ---------------------------------------------------------
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

hurr.N.enp=ddply(enp.hurr,c("WY.start"),summarise,N.val=N(Key))
hur.tck2=subset(hur.tck,Key%in%enp.hurr$Key)

hur.dat.enp=subset(hur.dat,Key%in%enp.hurr$Key&WY%in%seq(2005,2018,1))


## Accumulated Cyclone Energy (ACE)
max.wind=ddply(hur.dat.enp,c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
ENP.ACE=ddply(max.wind,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))
ENP.ACE=rbind(ENP.ACE,data.frame(WY=c(2008,2010,2012,2015),ACE=0))

max.wind.all=ddply(subset(hur.dat,WY%in%seq(2005,2018,1)),c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
ACE.all=ddply(max.wind.all,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))

plot(ACE~WY,ACE.all,ylim=c(0,350))
with(ENP.ACE,points(WY,ACE,pch=21,bg="red"))
abline(v=WY(c(wilma.landfall.FL,irma.landfall.FL)))


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

vars3=c("Station.ID","Region","Date.EST","TP","TN","TOC","DOC")
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
#grab.wq2=merge(grab.wq2,fce.wq.sites[,c("Station.ID","Region")],"Station.ID")

## Annual trend analysis
WQ.WY=ddply(grab.wq2,c("Station.ID","Region","datasource","EuDist.frac","WY"),summarise,mean.TP=mean(TP,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),N.TN=N(TN),mean.DOC=mean(DOC,na.rm=T),N.DOC=N(DOC))
WQ.WY$TP.scn=with(WQ.WY,ifelse(N.TP<4,NA,mean.TP))
WQ.WY$TN.scn=with(WQ.WY,ifelse(N.TN<4,NA,mean.TN))
WQ.WY$DOC.scn=with(WQ.WY,ifelse(N.DOC<4,NA,mean.DOC))
WQ.WY$dec.yr=with(WQ.WY,WY+0.328);#decimal date for Apirl 30th 
WQ.WY$WY.f=as.factor(WQ.WY$WY)
WQ.WY$WY.f2=paste0("WY",WQ.WY$WY)
WQ.WY$Hurr=as.factor(with(WQ.WY,ifelse(WY%in%c(2006,2018),1,0)))
WQ.WY=subset(WQ.WY,Station.ID%in%site.val)

## Time Series plot
unique(WQ.WY$Station.ID)
site.val=c(paste0("SRS",c("1d",2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))
WQ.WY$Station.ID=factor(WQ.WY$Station.ID,level=site.val)
#tiff(filename=paste0(plot.path,"FigX_WQPlots.tiff"),width=6.5,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2.75,0.5,1))
#layout(matrix(1:12,6,2,byrow=F))
layout(matrix(1:24,6,4,byrow=F))

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
  if(i==4){text(x=1999,y=1,"Total Phosphorus (\u03BCg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)

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
  if(i==4){text(x=1999.75,y=0.1,"Total Nitrogen (mg L\u207B\u00B9)",xpd=NA,srt=90,cex=1.75,adj=0)}
}
mtext(side=1,outer=T,"Water Year")
plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA)
legend.text=c("Annual Mean TP", "Annual Mean TN", "Grab Samples")
pt.cols=c("indianred1","dodgerblue1",adjustcolor("grey",0.75))
legend(0.5,0.75,legend=legend.text,pch =c(21,22,19),col=c("black","black",pt.cols[3]),pt.bg=pt.cols,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

##ANOVA (idea from from Davis et al 2018)
WQ.WY2=subset(WQ.WY,WY%in%seq(2005,2018,1))
#aov.TPcompare=aov(log(TP.scn)~Hurr+Station.ID+WY.f,subset(WQ.WY2,Region=="SRS"))
aov.TPcompare=aov(log(TP.scn)~Hurr+Station.ID+WY.f,WQ.WY2)
summary(aov.TPcompare)

#test assumptions
layout(matrix(1:4,2,2,byrow=F))
plot(aov.TPcompare)
dev.off()

# normality of residuals
shapiro.test(residuals(aov.TPcompare));
# homogenity of variances for each group
leveneTest(log(TP.scn)~Station.ID,WQ.WY); 
leveneTest(log(TP.scn)~WY.f,WQ.WY)
leveneTest(log(TP.scn)~Hurr,WQ.WY)

# Tukey HSD
tukey.station=TukeyHSD(aov.TPcompare,"Station.ID",ordered=T)
tukey.station.rslt=data.frame(tukey.station$Station.ID)
tukey.station.rslt$Station.ID=row.names(tukey.station.rslt)
tukey.station.rslt=cldList(p.adj~Station.ID,tukey.station.rslt,threshold = 0.05)
tukey.station.rslt=tukey.station.rslt[match(site.val,tukey.station.rslt$Group),]
tukey.station.rslt

x=boxplot(TP.scn~Station.ID,WQ.WY2,log="y",ylim=c(1,50),outline=F)
text(1:6,x$stats[5,]+5,toupper(tukey.station.rslt$Letter))

tukey.WY=TukeyHSD(aov.TPcompare,"WY.f",ordered=T)
tukey.WY.rslt=data.frame(tukey.WY$WY.f)
tukey.WY.rslt$WY=row.names(tukey.WY.rslt)
tukey.WY.rslt=cldList(p.adj~WY,tukey.WY.rslt,threshold = 0.05,remove.zero = F)
tukey.WY.rslt=tukey.WY.rslt[match(seq(2005,2018,1),tukey.WY.rslt$Group),]
tukey.WY.rslt

x=boxplot(TP.scn~WY,WQ.WY2,log="y",ylim=c(1,50),outline=F)
text(1:nrow(tukey.WY.rslt),x$stats[5,]+5,toupper(tukey.WY.rslt$Letter))

##
aov.TNcompare=aov(log(TN.scn)~Hurr+Station.ID+WY.f,WQ.WY2)
summary(aov.TNcompare)

#test assumptions
layout(matrix(1:4,2,2,byrow=F))
plot(aov.TNcompare)
dev.off()

# normality of residuals
shapiro.test(residuals(aov.TNcompare));
# homogenity of variances for each group
leveneTest(log(TN.scn)~Station.ID,WQ.WY); 
leveneTest(log(TN.scn)~WY.f,WQ.WY)
leveneTest(log(TN.scn)~Hurr,WQ.WY)

# Tukey HSD
tukey.station.TN=TukeyHSD(aov.TNcompare,"Station.ID",ordered=T)
tukey.station.TN.rslt=data.frame(tukey.station.TN$Station.ID)
tukey.station.TN.rslt$Station.ID=row.names(tukey.station.TN.rslt)
tukey.station.TN.rslt=cldList(p.adj~Station.ID,tukey.station.TN.rslt,threshold = 0.05)
tukey.station.TN.rslt=tukey.station.TN.rslt[match(site.val,tukey.station.TN.rslt$Group),]
tukey.station.TN.rslt

x=boxplot(TN.scn~Station.ID,WQ.WY2,log="y",ylim=c(0.1,5),outline=F)
text(1:11,x$stats[5,]+0.2,toupper(tukey.station.TN.rslt$Letter))

tukey.WY.TN=TukeyHSD(aov.TNcompare,"WY.f",ordered=T)
tukey.WY.TN.rslt=data.frame(tukey.WY.TN$WY.f)
tukey.WY.TN.rslt$WY=row.names(tukey.WY.TN.rslt)
tukey.WY.TN.rslt=cldList(p.adj~WY,tukey.WY.TN.rslt,threshold = 0.05,remove.zero = F)
tukey.WY.TN.rslt=tukey.WY.TN.rslt[match(seq(2005,2018,1),tukey.WY.TN.rslt$Group),]
tukey.WY.TN.rslt

x=boxplot(TN.scn~WY,WQ.WY2,log="y",ylim=c(0.1,5),outline=F)
text(1:nrow(tukey.WY.TN.rslt),x$stats[5,]+0.25,toupper(tukey.WY.TN.rslt$Letter))

## Pre-Post Storms (use autosample ... kruskal wallis test between pre-post like stage)

WY.storms=data.frame(WY=c(WY(c(wilma.landfall.FL,irma.landfall.FL)),WY(c(wilma.landfall.FL,irma.landfall.FL))-1),Period=c(rep("Post",2),rep("Pre",2)),Storm=rep(c("Wilma","Irma"),2))
#grab.wq3=merge(grab.wq2,analysis.periods,by.x="Date.EST",by.y="DATE")
grab.wq3=merge(grab.wq2,WY.storms,"WY")

plot(TP~Date.EST,subset(grab.wq3,Station.ID=="SRS1d"&WY%in%c(2017,2018,2019)))
abline(v=irma.landfall.FL)


#Climate CSF
WQ.WY2=merge(WQ.WY2,ENP.ACE,'WY')

WQ.WY2.TP=subset(WQ.WY2,is.na(TP.scn)==F)
subset(WQ.WY2.TP,is.na(Station.ID)==T)
##########BUILD MULTILEVEL MODELS: Simple 2-LEVEL MODEL------
#We want to predict annual TP
#Station.ID are the repeated unit 
#so you have to account for the non-independence of observations on the same Station

mnull<-lme(TP.scn~1,data=WQ.WY2.TP,random=~1|Station.ID,method="ML")

#How much does TP depend ACE?
#RI = Stations have different intercepts (biologically - they differ in average TP)
#RSI = Stations have different intercepts and slopes (biologically - they differ in the relationship between ACE and TP recorded across years)
mRI<-lme(TP.scn~ACE,data=WQ.WY2.TP,random=~1|Station.ID,method="ML")

#observations that occur in the same year may also be correlated, 
#so try random intercepts for year as well
mRI.2<-lme(TP.scn~ACE,data=WQ.WY2.TP,random=list(~1|Station.ID,~1|WY.f),method="ML")

#random slopes and intercepts model
mRSI<-lme(TP.scn~ACE,data=WQ.WY2.TP,random=~WY|Station.ID,method="ML")

#compare the three models: lowest AICc wins; difference <2 is a tie
mods.AICc=AICc(mnull, mRI, mRI.2, mRSI)
outer(mods.AICc$AICc,mods.AICc$AICc,"-")
mods.AICc
subset(mods.AICc,AICc==min(mods.AICc$AICc))

#null model is tied with mRI
# random intercept with station (not year)...makes sense since we have a strong gradient
# ACE explains some variation in TP relative to null model.

Anova(mRI,type=3)
# Chi-Square test for influence of ACE on TP fitted with type 3 error
# ACE does not significantly influence TP


#You could also get the p-value via a likelihood ratio test against the null model
#remember these models have to be nested - with the more complex model listed first
anova(mRI,mnull)

#estimated slope
summary(mRI)

rsquared(mRI)

#######MODEL SELECTION EXERCISE: NONLINEARITY ----------
#nonlinear Climate Sensitivity Functions (CSFs)
#Does TP vary nonlinearly with the hurricane intensity index (ACE)?

#Nonlinear responses to climate drivers are a way to assess
#whether TP will be sensitive (or not) to increasing VARIABILITY in Hurricane intensity
#Nonlinearities indicate sensitivity to change in variance of climate predictor
#Linear suggests no sensitivity to to change in variance of climate predictor

m1<-lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=~1|Station.ID,method="ML")
m2<-update(m1, .~. +I(ACE^2))
m3<-update(m2, .~. +I(ACE^3))
m3poly<-update(m3, .~poly(ACE,3))

nlmod.AICc=AICc(m1, m2, m3, m3poly)
nlmod.AICc
outer(nlmod.AICc$AICc,nlmod.AICc$AICc,"-")
subset(nlmod.AICc,AICc==min(nlmod.AICc$AICc))

#m2 is the best

Anova(m2,type=3)
#estimated slope
summary(m2)
rsquared(m2)

#######CHOOSING A VARIANCE-COVARIANCE MATRIX----------------

m1AR<-lme(log(TP.scn)~ACE,data=WQ.WY2.TP,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR<-update(m1AR, .~. +I(ACE^2))
m3AR<-update(m2AR, .~. +I(ACE^3))

#COMPARE MODELS based on AICc
nlmod.AICc2=AICc(m1,m1AR,m2,m2AR,m3,m3AR)
nlmod.AICc2

outer(nlmod.AICc2$AICc,nlmod.AICc2$AICc,"-")
subset(nlmod.AICc2,AICc==min(nlmod.AICc2$AICc))

anova(m2,m2AR)
summary(m2AR)
rsquared(m2AR)

## Check if SRS and TS models are different
m1AR.SRS<-lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="SRS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.SRS<-update(m1AR.SRS, .~. +I(ACE^2))
m3AR.SRS<-update(m2AR.SRS, .~. +I(ACE^3))

AICc(m1AR.SRS,m2AR.SRS,m3AR.SRS)
summary(m2AR.SRS)
rsquared(m2AR.SRS)
plot(m2AR.SRS)
qqnorm(residuals(m2AR.SRS));qqline(residuals(m2AR.SRS))

#
m1AR.TS<-lme(log(TP.scn)~ACE,data=subset(WQ.WY2.TP,Region=="TS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TS<-update(m1AR.TS, .~. +I(ACE^2))
m3AR.TS<-update(m2AR.TS, .~. +I(ACE^3))

AICc(m1AR.TS,m2AR.TS,m3AR.TS)
summary(m2AR.TS)
rsquared(m2AR.TS)
plot(m2AR.TS)
qqnorm(residuals(m2AR.TS));qqline(residuals(m2AR.TS))

par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2.75,0.5,1))
layout(matrix(1:2,2,1))
visreg(m2AR.SRS,"ACE",type="conditional",points.par=list(cex=1.2,col="black"))
visreg(m2AR.TS,"ACE",type="conditional",points.par=list(cex=1.2,col="dodgerblue1"))
dev.off()
###########EVALUATE MODEL ASSUMPTIONS---------------------

#1) Homogeneity of Variances (of best model)
plot(m2AR)

#2) Normality of Residuals (of best model)
qqnorm(residuals(m2AR));qqline(residuals(m2AR))

hist(residuals(m2AR))
shapiro.test(residuals(m2AR));#is significantly different from normal

######GET P-VALUES AND COEFFICIENT ESTIMATES WITH 95% CONFIDENCE INTERVALS-----------
#Likelihood Ratio X2 test
Anova(m2AR,type=3)

#F-tests (not as good as above for ML models, but here shows similar result)
anova.lme(m2AR,type = "marginal", adjustSigma = F)

intervals(m2AR, level = 0.95)

summary(m2AR)
rsquared(m2AR)

#sets graph size
visreg(m2AR,"ACE",type="conditional",points.par=list(cex=1.2,col="black"))
ylim.val=c(5,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,200);by.x=50;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

visreg(m2AR,"ACE",type="conditional",points.par=list(cex=1.2,col="black"),ylim=log(ylim.val),xlim=xlim.val,yaxt="n",ylab=NA,xaxt="n",xlab=NA)
axis_fun(2,log(ymaj),log(ymin),ymaj)
axis_fun(1,xmaj,xmin,xmaj)
mtext(side=1,line=2,"Accumulated Cyclone Energy (x10\u2074 knots\u00B2 Yr\u207B\u00B9)")
mtext(side=2,line=2,"Annual Mean TP (\u03BCg L\u207B\u00B9)")


vis.dat=visreg(m2AR,"ACE",type="conditional",plot=F)

plot(TP.scn~ACE,WQ.WY2.TP,log="y")
with(vis.dat$fit,lines(ACE,exp(visregFit)))
with(vis.dat$fit,lines(ACE,exp(visregLwr),col="red"))
with(vis.dat$fit,lines(ACE,exp(visregUpr),col="red"))


## 
WQ.WY2.TN=subset(WQ.WY2,is.na(TN.scn)==F)
subset(WQ.WY2.TN,is.na(Station.ID)==T)

#######MODEL SELECTION EXERCISE: NONLINEARITY ----------
#nonlinear Climate Sensitivity Functions (CSFs)
#Does TN vary nonlinearly with the hurricane intensity index (ACE)?

#Nonlinear responses to climate drivers are a way to assess
#whether TP will be sensitive (or not) to increasing VARIABILITY in Hurricane intensity
#Nonlinearities indicate sensitivity to change in variance of climate predictor
#Linear suggests no sensitivity to to change in variance of climate predictor

m1<-lme(log(TN.scn)~ACE,data=WQ.WY2.TN,random=~1|Station.ID,method="ML")
m2<-update(m1, .~. +I(ACE^2))
m3<-update(m2, .~. +I(ACE^3))

nlmod.AICc=AICc(m1, m2, m3)
nlmod.AICc
outer(nlmod.AICc$AICc,nlmod.AICc$AICc,"-")
subset(nlmod.AICc,AICc==min(nlmod.AICc$AICc))

#m1 is the best

Anova(m1,type=3)
#estimated slope
summary(m1)
rsquared(m1)

#1) Homogeneity of Variances (of best model)
plot(m1)

#2) Normality of Residuals (of best model)
qqnorm(residuals(m1));qqline(residuals(m1))

######CHOOSING A VARIANCE-COVARIANCE MATRIX----------------

m1AR<-lme(log(TN.scn)~ACE,data=WQ.WY2.TN,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR<-update(m1AR, .~. +I(ACE^2))
m3AR<-update(m2AR, .~. +I(ACE^3))

#COMPARE MODELS based on AICc
nlmod.AICc2=AICc(m1,m1AR,m2,m2AR,m3,m3AR)
nlmod.AICc2

outer(nlmod.AICc2$AICc,nlmod.AICc2$AICc,"-")
subset(nlmod.AICc2,AICc==min(nlmod.AICc2$AICc))

# m1AR is the best

anova(m1,m1AR)
summary(m1AR)
rsquared(m1AR)

## 
## Check if SRS and TS models are different
m1AR.SRS<-lme(log(TN.scn)~ACE,data=subset(WQ.WY2.TN,Region=="SRS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.SRS<-update(m1AR.SRS, .~. +I(ACE^2))
m3AR.SRS<-update(m2AR.SRS, .~. +I(ACE^3))

AICc(m1AR.SRS,m2AR.SRS,m3AR.SRS)
summary(m2AR.SRS)
rsquared(m2AR.SRS)
plot(m2AR.SRS)
qqnorm(residuals(m2AR.SRS));qqline(residuals(m2AR.SRS))
shapiro.test(residuals(m2AR.SRS))

#
m1AR.TS<-lme(log(TN.scn)~ACE,data=subset(WQ.WY2.TN,Region=="TS"),random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR.TS<-update(m1AR.TS, .~. +I(ACE^2))
m3AR.TS<-update(m2AR.TS, .~. +I(ACE^3))

AICc(m1AR.TS,m2AR.TS,m3AR.TS)
summary(m2AR.TS)
rsquared(m2AR.TS)
plot(m2AR.TS)
qqnorm(residuals(m2AR.TS));qqline(residuals(m2AR.TS))
shapiro.test(residuals(m2AR.TS))

par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2.75,0.5,1))
layout(matrix(1:2,2,1))
visreg(m2AR.SRS,"ACE",type="conditional",points.par=list(cex=1.2,col="black"))
visreg(m2AR.TS,"ACE",type="conditional",points.par=list(cex=1.2,col="dodgerblue1"))
dev.off()

###########EVALUATE MODEL ASSUMPTIONS---------------------

#1) Homogeneity of Variances (of best model)
plot(m1AR)

#2) Normality of Residuals (of best model)
qqnorm(residuals(m1AR));qqline(residuals(m1AR))

hist(residuals(m1AR))
shapiro.test(residuals(m1AR));#is significantly different from normal

######GET P-VALUES AND COEFFICIENT ESTIMATES WITH 95% CONFIDENCE INTERVALS-----------
#Likelihood Ratio X2 test
Anova(m1AR,type=3)

#F-tests (not as good as above for ML models, but here shows similar result)
anova.lme(m1AR,type = "marginal", adjustSigma = F)

intervals(m1AR, level = 0.95)

summary(m1AR)
rsquared(m1AR)

#sets graph size
visreg(m1AR,"ACE",type="conditional",points.par=list(cex=1.2,col="black"))

## 
WQ.WY2.DOC=subset(WQ.WY2,is.na(DOC.scn)==F)
subset(WQ.WY2.DOC,is.na(Station.ID)==T)

m1<-lme(log(DOC.scn)~ACE,data=WQ.WY2.DOC,random=~1|Station.ID,method="ML")
m2<-update(m1, .~. +I(ACE^2))
m3<-update(m2, .~. +I(ACE^3))

m1AR<-lme(log(DOC.scn)~ACE,data=WQ.WY2.DOC,random=~1|Station.ID,correlation=corAR1(form=~WY),method="ML")
m2AR<-update(m1AR, .~. +I(ACE^2))
m3AR<-update(m2AR, .~. +I(ACE^3))

#COMPARE MODELS based on AICc
nlmod.AICc2=AICc(m1,m1AR,m2,m2AR,m3,m3AR)
nlmod.AICc2

outer(nlmod.AICc2$AICc,nlmod.AICc2$AICc,"-")
subset(nlmod.AICc2,AICc==min(nlmod.AICc2$AICc))

# m1 is the best

anova(m1,m1AR)
summary(m1AR)
rsquared(m1AR)

#1) Homogeneity of Variances (of best model)
plot(m1AR)

#2) Normality of Residuals (of best model)
qqnorm(residuals(m1AR));qqline(residuals(m1AR))

shapiro.test(residuals(m1AR))

#Likelihood Ratio X2 test
Anova(m1AR,type=3)

visreg(m1AR,"ACE",type="conditional",points.par=list(cex=1.2,col="black"))
