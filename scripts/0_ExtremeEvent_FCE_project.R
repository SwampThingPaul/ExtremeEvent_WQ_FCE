
## Frontiers in Science 
## (Effects of Extreme Weather Events on Coastal Carbon and Nutrient Cycling)
## 
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(plyr)
library(reshape)
library(zoo);
library(lubridate);

library(dunn.test)
library(rcompanion)
library(smatr)
library(hydrostats)
library(tidyr)
library(gvlma)
library(mblm)
library(trend)
library(segmented)
library(vegan)

#GIS libraries
library(tmap)
library(tmaptools)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)


#Custom Functions
source("hydrostats_pj.r")

source("D:/CommonlyUsedFunctions.r")
#the commonly used functions is also here
#source("https://github.com/SwampThingPaul/analyst_helper/blob/3bc5585bcbe2af65fdf36172dd4047e4791f594d/CommonlyUsedFunctions.r")
RBFlash <- function (m, na.rm=TRUE) {
  # from https://rdrr.io/github/leppott/ContDataQC/src/R/RBIcalc.R
  sum(abs(diff(m)), na.rm=na.rm)/sum(m, na.rm=na.rm)
}
#Paths
setwd("D:/UF/LTER/Projects/ExtremeEvent_WQ_FCE")

paths=paste0(getwd(),c("/Exports/","/Plots/","/Data/","/GIS"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]


GIS.path.gen="D:/_GISData"
#shell.exec(system.file(getwd()))

#Helper variables 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


# ------------------------------------------------------------------------
enp=readOGR(paste0(GIS.path.gen,"/SFER_GIS_Geodatabase.gdb"),"ENP")
fce.sites=spTransform(readOGR(GIS.path,"ltersites_utm"),utm17)

enp.shoreline=spTransform(readOGR(GIS.path,"ENP_shoreline"),utm17)
enp.slough=spTransform(readOGR(GIS.path,"sloughs_utm"),utm17)
enp.slough.clip=gIntersection(enp.shoreline,enp.slough)
sfwmd.mon=readOGR(paste0(GIS.path.gen,"/SFWMD_Monitoring_20180829"),"Environmental_Monitoring_Stations")
sfwmd.mon=spTransform(sfwmd.mon,utm17)
sfwmd.mon$UTMx=coordinates(sfwmd.mon)[,1]
sfwmd.mon$UTMy=coordinates(sfwmd.mon)[,2]
#sfwmd.mon=subset(sfwmd.mon,ACTIVITY_S%in%c("Stage","Well")&STATUS=="Active")
sfwmd.mon=subset(sfwmd.mon,STATUS=="Active")
unique(sfwmd.mon$ACTIVITY_S)

sfwmd.mon.wq=subset(sfwmd.mon,ACTIVITY_S=="Surface Water Grab")
sfwmd.mon.wx=subset(sfwmd.mon,ACTIVITY_S%in%c("Weather","Rain"))
sfwmd.mon.q=subset(sfwmd.mon,ACTIVITY_S%in%c("Flow")&STATUS=="Active")
# Storm Data --------------------------------------------------------------
pre.post.dur=2
wilma.landfall.FL=date.fun("2005-10-24")
wilma.period=date.fun(c(wilma.landfall.FL-duration(pre.post.dur,"months"),wilma.landfall.FL+duration(pre.post.dur,"months")))#date.fun(c("2003-10-01","2005-01-01"))
irma.landfall.FL=date.fun("2017-09-10")
irma.period=date.fun(c(irma.landfall.FL-duration(pre.post.dur,"months"),irma.landfall.FL+duration(pre.post.dur,"months")))#date.fun(c("2016-10-01","2018-01-01"))
analysis.period.pre.wilma=date.fun(seq(wilma.period[1],wilma.landfall.FL,"1 days"))
analysis.period.post.wilma=date.fun(seq(wilma.landfall.FL+ddays(1),wilma.period[2],"1 days"))
analysis.period.pre.irma=date.fun(seq(irma.period[1],irma.landfall.FL,"1 days"))
analysis.period.post.irma=date.fun(seq(irma.landfall.FL+ddays(1),irma.period[2],"1 days"))
analysis.periods=data.frame(DATE=c(analysis.period.pre.wilma,analysis.period.post.wilma,analysis.period.pre.irma,analysis.period.post.irma),
                            Storm=c(rep("Wilma",length(c(analysis.period.pre.wilma,analysis.period.post.wilma))),rep("Irma",length(c(analysis.period.pre.irma,analysis.period.post.irma)))),
                            Period=c(rep("Pre",length(analysis.period.pre.wilma)),rep("Post",length(analysis.period.post.wilma)),rep("Pre",length(analysis.period.pre.irma)),rep("Post",length(analysis.period.post.irma))))


# Rainfall  ---------------------------------------------------------------
#source(paste0(getwd(),"/scripts/1_ExtremeEvent_Rainfall.r"),echo=T)


# Water Level -------------------------------------------------------------
#source(paste0(getwd(),"/scripts/2_ExtremeEvent_stage.r"),echo=T)


# Weather Data ------------------------------------------------------------
#source(paste0(getwd(),"/scripts/3_ExtremeEvent_weather.r"),echo=T)

# Discharge & Load --------------------------------------------------------
#source(paste0(getwd(),"/scripts/4_ExtremeEvent_flowload.r"),echo=T)

# Water Quality Gradient --------------------------------------------------
tmap_mode("view")
tm_shape(enp.slough.clip)+tm_polygons(alpha=0.5)+
  tm_shape(fce.sites)+tm_dots(col="yellow")+
  tm_shape(sfwmd.mon.wq)+tm_dots(col="red")


sdate=as.Date("1999-05-01");# Begining of Florida WY2005
edate=as.Date("2018-04-30");# End of Florida WY2018
wq.params=data.frame(param=c("TP","TN","NOx","TKN","DOC","TOC"),Test.Number=c(25,80,18,21,89,100))

slough.wq.sites=data.frame(WQSite=c("NE1","P33","P36","P35","TSB","P37"),
                    Region=c(rep("SRS",4),rep("TS",2)));#FLAB sites from WMD spreadsheet
slough.wq=SFWMD.DBHYDRO.Data.WQ(sdate,edate,slough.wq.sites$WQSite,wq.params$Test.Number)
slough.wq=subset(slough.wq,Collection.Method=="G")
slough.wq$AbsMDL=with(slough.wq,ifelse(Value<0,abs(Value),Value))
slough.wq=merge(slough.wq,wq.params,"Test.Number")
slough.wq=merge(slough.wq,slough.wq.sites,by.x="Station.ID",by.y="WQSite")
slough.wq.xtab=cast(slough.wq,Station.ID+Region+Date.EST~param,value="AbsMDL",mean)
slough.wq.xtab$TN.final=with(slough.wq.xtab,SFWMD.TN.Combine(NOx,TKN,TN))
slough.wq.xtab$TOC=NA

vars=c("Station.ID","Region","Date.EST","TP","TN.final","TOC")
slough.wq.xtab=slough.wq.xtab[,vars]
slough.wq.xtab=rename(slough.wq.xtab,c("TN.final"="TN"))


#coast.wq.sites2=data.frame(WQSite=c("FLAB33","FLAB34","FLAB35","FLAB38","FLAB37","FLAB36","FLAB39","FLAB40","FLAB11"),
#                           Region2=c(rep("SRS",8),"TS"))
coast.wq.sites2=data.frame(WQSite=c("FLAB38","FLAB39","FLAB40","FLAB11"),
                           Region2=c(rep("SRS",3),"TS"))
coast.wq.xtab3=read.csv(paste0(data.path,"SFWMD/CWQMN_thru_April2018_For_Analysis.csv"))
coast.wq.xtab3=merge(coast.wq.xtab3,coast.wq.sites2,by.x="Station.ID",by.y="WQSite")
coast.wq.xtab3$Date.EST=date.fun(coast.wq.xtab3$Collection.Date,form="%m/%d/%Y")
coast.wq.xtab3$TP=(coast.wq.xtab3$CTP.P.um.l*P.mw)*0.001
coast.wq.xtab3$TN=coast.wq.xtab3$TN.mg.l
coast.wq.xtab3$TOC=(coast.wq.xtab3$TOC.um.l*C.mw)*0.001

vars2=c("Station.ID","Region2","Date.EST","TP","TN","TOC")
coast.wq.xtab3=coast.wq.xtab3[,vars2]
coast.wq.xtab3=rename(coast.wq.xtab3,c("Region2"="Region"))

wmd.wq=rbind(slough.wq.xtab,coast.wq.xtab3)
wmd.wq$DOC=NA
head(wmd.wq)
## FCE data (grab data)
fce.grahl=read.csv(paste0(data.path,"/WQ/LT_ND_Grahl_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.losada=read.csv(paste0(data.path,"/WQ/LT_ND_Losada_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rondeau=read.csv(paste0(data.path,"/WQ/LT_ND_Rondeau_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rubio=read.csv(paste0(data.path,"/WQ/LT_ND_Rubio_002.txt"),na.strings=c("-9999","-9999.00","-9999.000"))

names(fce.grahl)
names(fce.losada)
names(fce.rondeau)
names(fce.rubio)
fce.grahl=rename(fce.grahl,c("NandN"="N.N"))
fce.wq=rbind(fce.grahl,fce.losada,fce.rondeau,fce.rubio)
fce.wq$Date.EST=date.fun(as.character(fce.wq$DATE),form="%Y-%m-%d")
#subset(fce.wq,is.na(Date.EST)==T)

unique(fce.wq$SITENAME)
#fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,"6b","6a","7a","7b"))),
#                     Station.ID=c(c(rep("SRS1",3),paste0("SRS",2:6)),c("TS/PH1","TS/PH1","TS/PH2","TS/PH6a","TS/PH6b","TS/PH7a","TS/PH7b")),
#                     Region=c(rep("SRS",8),rep("TS",7)))
fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                     Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                     Region=c(rep("SRS",8),rep("TS",8)))

fce.wq=merge(fce.wq,fce.wq.sites,"SITENAME")
unique(fce.wq$SITENAME)

fce.wq$TP=round((fce.wq$TP*P.mw)*0.001,4)
fce.wq$TN=round((fce.wq$TN*N.mw)*0.001,2)
fce.wq$TOC=round((fce.wq$TOC*C.mw)*0.001,2)
fce.wq$DOC=round((fce.wq$DOC*C.mw)*0.001,2)

vars3=c("Station.ID","Region","Date.EST","TP","TN","TOC","DOC")
fce.wq=fce.wq[,vars3]
head(fce.wq)

grab.wq=rbind(fce.wq,wmd.wq)

### 
coast.wq.sites2=rename(coast.wq.sites2,c("Region2"="Region"))
wmd.sites.all=rbind(slough.wq.sites,coast.wq.sites2)

# spatial data
wmd.wq.sp=merge(subset(sfwmd.mon.wq,STATION%in%wmd.sites.all$WQSite),wmd.sites.all,by.x="STATION",by.y="WQSite")
wmd.wq.sp$datasource="SFWMD"
unique(wmd.wq.sp$Region)

fce.sp.site.list=c(paste0("SRS-",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/Ph-",c("1a","1b","2a","2b",2,3,"6b","6a","7a","7b")))
fce.sites2=subset(fce.sites,SITE%in%fce.sp.site.list)
fce.sites2=merge(fce.sites2,data.frame(SITE=fce.sp.site.list,STATION=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b","2a","2b",2,3,"6b","6a","7a","7b"))),Region=c(rep("SRS",8),rep("TS",10))),"SITE")
fce.sites2$datasource="FCE"
wq.sites.sp=rbind(wmd.wq.sp[,c("STATION","Region","datasource")],fce.sites2[,c("STATION","Region","datasource")])

TS2.locdat=data.frame(STATION="TS/PH2",Region="TS",datasource="FCE",LONG=-80.607,LAT=25.404)
TS2=spTransform(SpatialPointsDataFrame(coords=TS2.locdat[,c("LONG","LAT")],data=TS2.locdat[,c("STATION","Region","datasource")],proj4string=CRS("+init=epsg:4326")),utm17)
wq.sites.sp=rbind(wq.sites.sp,TS2)

tm_shape(wmd.wq.sp)+tm_dots(col="red")
tm_shape(fce.sites)+tm_dots()
tm_shape(wq.sites.sp)+tm_dots(col="datasource",size=0.1,palette=c("blue","green"))

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
head(wq.sites.sp@data)

grab.wq2=merge(grab.wq,wq.sites.sp@data[,c("STATION","datasource","EuDist.frac")],by.x="Station.ID",by.y="STATION",all.x=T)
head(grab.wq2)
grab.wq2$month=format(grab.wq2$Date.EST,"%m")
grab.wq2$CY=as.numeric(format(grab.wq2$Date.EST,"%Y"))
grab.wq2$WY=WY(grab.wq2$Date.EST)
#grab.wq2=subset(grab.wq2,WY%in%seq(2005,2018,1))

## Annual trend analysis
WQ.WY=ddply(grab.wq2,c("Station.ID","datasource","EuDist.frac","WY"),summarise,mean.TP=mean(TP,na.rm=T),N.TP=N(TP))
WQ.WY$TP.scn=with(WQ.WY,ifelse(N.TP<4,NA,mean.TP))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS1d"))

with(subset(WQ.WY,Station.ID=="SRS1d"),cor.test(TP.scn,WY,method="kendall"))
pettitt.test(subset(WQ.WY,Station.ID=="SRS1d")$TP.scn)
mod=lm(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS1d"))
seg.mod=segmented(mod,seg.Z=~WY)
summary(seg.mod)
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS1d"))
plot(seg.mod,add=T)

with(subset(WQ.WY,Station.ID=="SRS2"),cor.test(TP.scn,WY,method="kendall"))
pettitt.test(subset(WQ.WY,Station.ID=="SRS2")$TP.scn)
mod=lm(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS2"&is.na(TP.scn)==F))
#seg.mod=segmented(mod,seg.Z=~WY)
#summary(seg.mod)
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS2"))
#plot(seg.mod,add=T)

with(subset(WQ.WY,Station.ID=="SRS3"),cor.test(TP.scn,WY,method="kendall"))
pettitt.test(subset(WQ.WY,Station.ID=="SRS3")$TP.scn)
mod=lm(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS3"&is.na(TP.scn)==F))
#seg.mod=segmented(mod,seg.Z=~WY)
#summary(seg.mod)
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS3"))
#plot(seg.mod,add=T)

plot(TP.scn~WY,subset(WQ.WY,Station.ID%in%paste0("SRS",c("1a","1c","1b","1d"))))
with(subset(WQ.WY,Station.ID%in%paste0("SRS",c("1a","1c","1b","1d"))),text(x=WY,y=TP.scn+0.001,label=substring(Station.ID,4,6)))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS1a"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS1d"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS2"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS3"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS4"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS5"))
plot(TP.scn~WY,subset(WQ.WY,Station.ID=="SRS6"))


summary(lm(TP.scn~WY*EuDist.frac,WQ.WY))



## FCE data
grab.wq2.fce=subset(grab.wq2,datasource=="FCE")
unique(grab.wq2.fce$Station.ID)

vars=c("Date","WQSite","Data.Value.m","depth.cm","hydroperiod.90dfreq")
wq.wl=merge(grab.wq2.fce,WL.dat[,vars],by.x=c("Station.ID","Date.EST"),by.y=c("WQSite","Date"))
wq.wl=merge(wq.wl,WL.dat.xtab.melt,by.x=c("Date.EST","Region"),by.y=c("Date","Region"))
wq.wl=merge(wq.wl,rf.dat2[,c("Date.EST","TRF.m","cum.rf")],"Date.EST")

plot(TP~cum.rf,subset(wq.wl,Region=="SRS"))
names(wq.wl)

## PCA Analysis
my.rda=rda(tmp.dat[,-1:-4],scale=T)
biplot(my.rda)

eig <- my.rda$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)

#tiff(filename=paste0(plot.path,"/PCA_screplots.tiff"),width=6.75,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,1.5,1,1),oma=c(3.5,2,0.25,0.25),mgp=c(3,1,0));
layout(matrix(seq(1,2,1),1,2,byrow=T));

ylim.val=c(0,105);by.y=20;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)
points(x,eig.pca.all$cumvariance,pch=21,bg="indianred1",cex=1.25)
axis_fun(1,line=-0.7,x,x,seq(1,9,1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=2,"Percentage of Variances")
leg.x=x[1,]+(x[length(vars)-1,]-x[1,])/2
leg.y=ylim.val[1]-25
legend.text=c("Absolute","Cumulative")
pt.col=c("grey","indianred1")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1.55,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

ylim.val=c(0,3);by.y=1;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n",add=T)
axis_fun(1,line=-0.7,x,x,seq(1,9,1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=1.5,"Eigenvalue")
dev.off()

plot(c(-2,2),c(-2,2),type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
x=ordiellipse(my.rda,group=tmp.dat$Region,draw="polygon",label=F,col=c("dodgerblue1","indianred1"),border=T,cex=0.8)

labs.rownames=merge(data.frame(vars=rownames(var.coord)),rownames.xwalk,sort=F)$short
scrs=scores(my.rda,display=c("sites","species"))

plot(c(-2,2),c(-2,2),type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
x=ordiellipse(my.rda,group=tmp.dat$Region,draw="polygon",label=F,col=c("indianred1","dodgerblue1"),border=T,cex=0.8)

labs.val=c("SRS","TS")
ellipse.center=data.frame()
for(i in 1:length(labs.val)){
  tmp=data.frame(Region=labs.val[i],pc1=as.numeric(x[[labs.val[i]]]$center[1]),pc2=as.numeric(x[[labs.val[i]]]$center[2]))
  ellipse.center=rbind(ellipse.center,tmp)
}
label_line=function(cen.x,cen.y,x2,y2,label,lty=1,cex=1,col="black",lwd=1,font=1,pos=NULL,pt.plot=T){
  if(pt.plot==F){NA}else{points(cen.x,cen.y,pch=19,col=col,cex=cex)}
  lines(c(cen.x,x2-(x2*0.06)),c(cen.y,y2-(y2*0.06)),lty=lty,col=col,lwd=lwd)
  text(x2,y2,label=label,cex=cex,col=col,font=font,pos=pos)
}

#tiff(filename=paste0(plot.path,"PCA_FCE.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.5,1.75),oma=c(2,1.75,0.25,0),mgp=c(3,1,0));
#layout(matrix(seq(1,2,1),1,2,byrow=T));

xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-3,3);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
with(tmp.dat,points(scrs$sites,pch=ifelse(tmp.dat$Region=="SRS",21,22),bg=adjustcolor(ifelse(tmp.dat$Region=="SRS","grey","white"),0.5),col=adjustcolor(ifelse(tmp.dat$Region=="SRS","grey30","black"),0.5),cex=1,lwd=0.5))
arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.1, angle = 15, code = 2,col="indianred1",lwd=1.5)
with(scrs,text(species[,1]+0.25,species[,2],labels=labs.rownames,cex=0.75))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"))
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"))
legend("topleft",legend=c("SRS","TS"),pch=c(21,22),lty=0,pt.bg=adjustcolor(c("grey","white"),0.5),col=adjustcolor(c("grey30","black"),0.5),pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()






###
grab.wq2.month=ddply(subset(grab.wq2,Station.ID!="P36"),c("Station.ID","Region","datasource","EuDist.frac","month","CY"),summarise,mean.TP=mean(TP*1000,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),mean.DOC=mean(DOC,na.rm=T))
#grab.wq2.month=ddply(subset(grab.wq2,datasource=="FCE"),c("Station.ID","Region","datasource","EuDist.frac","month","CY"),summarise,mean.TP=mean(TP*1000,na.rm=T))
#grab.wq2.month=ddply(grab.wq2,c("Station.ID","Region","datasource","EuDist.frac","month","CY"),summarise,mean.TP=mean(TP*1000,na.rm=T))
exclude.sites=c("SRS1a","SRS1c","P36","P35","FLAB39","FLAB38")
grab.wq2.month=ddply(subset(grab.wq2,!(Station.ID%in%exclude.sites)),c("Station.ID","Region","datasource","EuDist.frac","month","CY"),summarise,mean.TP=mean(TP*1000,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),mean.DOC=mean(DOC,na.rm=T))
grab.wq2.month$EuDist.frac=round(grab.wq2.month$EuDist.frac,2)
grab.wq2.month$mon.cy.date=date.fun(with(grab.wq2.month,paste(CY,month,"01",sep="-")))
grab.wq2.month$dec.year=decimal_date(grab.wq2.month$mon.cy.date)
grab.wq2.month=subset(grab.wq2.month,is.na(mean.TP)==F)
###

sites.val=ddply(subset(grab.wq2.month,Region=="SRS"&!(Station.ID%in%exclude.sites)),c("EuDist.frac","Station.ID"),summarise,N.val=N(mean.TP))#unique(subset(grab.wq2.month,Region=="SRS")$Station.ID)
par(family="serif",mar=c(1.5,3,0.5,1.5),oma=c(2,2,1,1));
layout(matrix(1:14,ncol=2))

for(i in 1:nrow(sites.val)){
  plot(mean.TP~mon.cy.date,subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),ylim=c(2,100),log="y")
  mtext(side=3,line=-1,sites.val$Station.ID[i],col="red",cex=0.8)
}




dates.val=ddply(grab.wq2.month,c("Region","mon.cy.date"),summarise,N.val=N(mean.TP))#sort(unique(grab.wq2.month$mon.cy.date))
dates.val=subset(dates.val,N.val>4)
N.cols=nrow(subset(dates.val,Region=="SRS"))
cols=grey.colors(N.cols)

plot(mean.TP~EuDist.frac,grab.wq2.month,type="n")
for(i in 1:N.cols){
  with(subset(grab.wq2.month,mon.cy.date==subset(dates.val,Region=="SRS")$mon.cy.date[i]),lines(EuDist.frac,mean.TP,col=adjustcolor(cols[i],0.5),lwd=2.5))
}



grab.wq2.WY=ddply(subset(grab.wq2,!(Station.ID%in%exclude.sites)),c("EuDist.frac","Station.ID","Region","datasource","WY"),summarise,mean.TP=mean(TP*1000,na.rm=T),N.TP=N(TP),mean.TN=mean(TN,na.rm=T),mean.DOC=mean(DOC,na.rm=T))
grab.wq2.WY.TP=subset(grab.wq2.WY,N.TP>4)
WYS=seq(2005,2018,1)
cols=grey.colors(length(WYS))

#tiff(filename=paste0(plot.path,"/SRS_TP_gradient.tiff"),width=10,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"/SRS_TP_gradient.png"),width=10,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,1,0.5,0.5),oma=c(2,2.5,1,1));
layout(matrix(c(1:9,10:11,rep(12,5),13:14),2,9,byrow=T))

xlim.val=date.fun(c("2004-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(2,200);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:nrow(sites.val)){
  plot(mean.TP~mon.cy.date,subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),ylim=ylim.val,xlim=xlim.val,log="y",yaxt="n",xaxt="n",type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),points(mon.cy.date,mean.TP,pch=21,cex=1,col=adjustcolor("black",0.5),bg=adjustcolor("grey",0.5)))
  mtext(side=3,line=-1,sites.val$Station.ID[i],cex=0.8)
  if(i==1){axis_fun(2,ymaj,ymin,ymaj,0.8)}else{axis_fun(2,ymaj,ymin,NA)}
  if(i==1){mtext(side=2,line=2,"TP (\u03BCg L\u207B\u00B9)")}
  axis_fun(1,line=-0.70,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.8)
  k=with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),lowess(mon.cy.date,mean.TP))
  lines(date.fun(as.POSIXct(k$x,origin="1970-01-01")),k$y,col=adjustcolor("red",0.5),lwd=2.5)
  }
plot(1:0,axes=F,type="n")
plot(1:0,axes=F,type="n")

ylim.val=c(0,50);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~EuDist.frac,grab.wq2.WY,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYS)){
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),lines(EuDist.frac,mean.TP,col=adjustcolor(cols[i],0.5),lwd=2.5))
  k=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),lowess(EuDist.frac,mean.TP,f=0.5))
  #lines(k$x,k$y,col=adjustcolor(cols[i],0.5),lwd=5)
  #ss=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),smooth.spline(EuDist.frac,mean.TP,df=7))
  #lines(ss)
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),points(EuDist.frac,mean.TP,pch=21,bg=adjustcolor(cols[i],0.75),cex=1.25))
}
axis_fun(1,xmaj,xmin,format(xmaj),1);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=2,line=2,"Annual Mean TP (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=2,"Fraction Distance Downstream")

legend_image <- as.raster(rev(adjustcolor(cols,0.5)), ncol=1)
plot(c(0,2),c(0,1),type = 'n', axes=F,xlab=NA, ylab=NA)
text(1,1.05,"WY",xpd=NA)
text(x=1.5, y = c(0,0.5,1), labels = c(2005,2012,2018))
rasterImage(legend_image, 0,0,1,1)
plot(1:0,axes=F,type="n")
dev.off()

sites.val=ddply(subset(grab.wq2.month,Region=="SRS"&!(Station.ID%in%exclude.sites)),c("EuDist.frac","Station.ID"),summarise,N.val=N(mean.TN))
#tiff(filename=paste0(plot.path,"/SRS_TN_gradient.tiff"),width=10,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,1,0.5,0.5),oma=c(2,2.5,1,1));
layout(matrix(c(1:9,10:11,rep(12,5),13:14),2,9,byrow=T))

xlim.val=date.fun(c("2004-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(0.1,10);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:nrow(sites.val)){
  plot(mean.TP~mon.cy.date,subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),ylim=ylim.val,xlim=xlim.val,log="y",yaxt="n",xaxt="n",type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),points(mon.cy.date,mean.TN,pch=21,cex=1,col=adjustcolor("black",0.5),bg=adjustcolor("grey",0.5)))
  mtext(side=3,line=-1,sites.val$Station.ID[i],cex=0.8)
  if(i==1){axis_fun(2,ymaj,ymin,ymaj,0.8)}else{axis_fun(2,ymaj,ymin,NA)}
  if(i==1){mtext(side=2,line=2,"TN (mg L\u207B\u00B9)")}
  axis_fun(1,line=-0.70,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.8)
  #k=with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),lowess(mon.cy.date,mean.TN,f=0.001))
  #lines(date.fun(as.POSIXct(k$x,origin="1970-01-01")),k$y,col=adjustcolor("red",0.5),lwd=2.5)
}
plot(1:0,axes=F,type="n")
plot(1:0,axes=F,type="n")

ylim.val=c(0,2);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~EuDist.frac,grab.wq2.WY,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYS)){
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),lines(EuDist.frac,mean.TN,col=adjustcolor(cols[i],0.5),lwd=2.5))
  k=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),lowess(EuDist.frac,mean.TP,f=0.1))
  #lines(k$x,k$y,col=adjustcolor(cols[i],0.5),lwd=5)
  #ss=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),smooth.spline(EuDist.frac,mean.TP,df=7))
  #lines(ss)
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),points(EuDist.frac,mean.TN,pch=21,bg=adjustcolor(cols[i],0.75),cex=1.25))
}
axis_fun(1,xmaj,xmin,format(xmaj),1);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=2,line=2,"Annual Mean TN (mg L\u207B\u00B9)")
mtext(side=1,line=2,"Fraction Distance Downstream")

legend_image <- as.raster(rev(adjustcolor(cols,0.5)), ncol=1)
plot(c(0,2),c(0,1),type = 'n', axes=F,xlab=NA, ylab=NA)
text(1,1.05,"WY",xpd=NA)
text(x=1.5, y = c(0,0.5,1), labels = c(2005,2012,2018))
rasterImage(legend_image, 0,0,1,1)
plot(1:0,axes=F,type="n")
dev.off()

sites.val=ddply(subset(grab.wq2.month,Region=="SRS"&!(Station.ID%in%exclude.sites)),c("EuDist.frac","Station.ID"),summarise,N.val=N(mean.DOC))
#tiff(filename=paste0(plot.path,"/SRS_DOC_gradient.tiff"),width=10,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(2,1,0.5,0.5),oma=c(2,2.5,1,1));
layout(matrix(c(1:9,10:11,rep(12,5),13:14),2,9,byrow=T))

xlim.val=date.fun(c("2004-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(1,100);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
for(i in 1:nrow(sites.val)){
  plot(mean.TP~mon.cy.date,subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),ylim=ylim.val,xlim=xlim.val,log="y",yaxt="n",xaxt="n",type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),points(mon.cy.date,mean.DOC,pch=21,cex=1,col=adjustcolor("black",0.5),bg=adjustcolor("grey",0.5)))
  mtext(side=3,line=-1,sites.val$Station.ID[i],cex=0.8)
  if(i==1){axis_fun(2,ymaj,ymin,ymaj,0.8)}else{axis_fun(2,ymaj,ymin,NA)}
  if(i==1){mtext(side=2,line=2,"DOC (mg L\u207B\u00B9)")}
  axis_fun(1,line=-0.70,xmaj,xmin,format(xmaj,"%m-%Y"),cex=0.8)
  k=with(subset(grab.wq2.month,Station.ID==sites.val$Station.ID[i]),lowess(mon.cy.date,mean.DOC,f=0.4))
  lines(date.fun(as.POSIXct(k$x,origin="1970-01-01")),k$y,col=adjustcolor("red",0.5),lwd=2.5)
}
plot(1:0,axes=F,type="n")
plot(1:0,axes=F,type="n")

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(mean.TP~EuDist.frac,grab.wq2.WY,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:length(WYS)){
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"&is.na(mean.DOC)==F),lines(EuDist.frac,mean.DOC,col=adjustcolor(cols[i],0.5),lwd=2.5,lty=2))
  k=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),lowess(EuDist.frac,mean.DOC,f=0.1))
  #lines(k$x,k$y,col=adjustcolor(cols[i],0.5),lwd=5)
  #ss=with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"),smooth.spline(EuDist.frac,mean.TP,df=7))
  #lines(ss)
  with(subset(grab.wq2.WY.TP,WY==WYS[i]&Region=="SRS"&is.na(mean.DOC)==F),points(EuDist.frac,mean.DOC,pch=21,bg=adjustcolor(cols[i],0.75),cex=1.25))
}
axis_fun(1,xmaj,xmin,format(xmaj),1);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=2,line=2,"Annual Mean DOC (mg L\u207B\u00B9)")
mtext(side=1,line=2,"Fraction Distance Downstream")

legend_image <- as.raster(rev(adjustcolor(cols,0.5)), ncol=1)
plot(c(0,2),c(0,1),type = 'n', axes=F,xlab=NA, ylab=NA)
text(1,1.05,"WY",xpd=NA)
text(x=1.5, y = c(0,0.5,1), labels = c(2005,2012,2018))
rasterImage(legend_image, 0,0,1,1)
plot(1:0,axes=F,type="n")
dev.off()


####
dates.val
grab.wq2.month.TP=ddply(subset(grab.wq2,!(Station.ID%in%exclude.sites)),c("Station.ID","Region","datasource","EuDist.frac","month","CY"),summarise,mean.TP=mean(TP*1000,na.rm=T),N.TP=N(TP),log.TP=log(mean.TP))
grab.wq2.month.TP$mon.cy.date=date.fun(with(grab.wq2.month.TP,paste(CY,month,"01",sep="-")))
grab.wq2.month.TP$dec.year=decimal_date(grab.wq2.month.TP$mon.cy.date)
grab.wq2.month.TP=subset(grab.wq2.month.TP,is.na(mean.TP)==F)
subset(grab.wq2.month.TP,N.TP>1)

#SRS
dates.val.srs=subset(dates.val,Region=="SRS")

SRS.tp.gradient=data.frame()
for(i in 1:nrow(dates.val.srs)){
  tmp.dat=subset(grab.wq2.month.TP,Region=="SRS"&mon.cy.date==dates.val.srs$mon.cy.date[i])
  #plot(mean.TP~EuDist.frac,tmp.dat,log="y")
  mblm.mod=mblm(log.TP~EuDist.frac,tmp.dat)
  
  tmp.rslt=data.frame(mon.cy.date=dates.val.srs$mon.cy.date[i],Region="SRS",N.sites=nrow(tmp.dat),slope.val=as.numeric(coefficients(mblm.mod)[2]),r2=summary.lm(mblm.mod)$r.squared,RSE=summary(mblm.mod)$sigma)
  SRS.tp.gradient=rbind(SRS.tp.gradient,tmp.rslt)
  print(i)
}
plot(slope.val~mon.cy.date,SRS.tp.gradient,type="l")
abline(v=c(irma.landfall.FL,wilma.landfall.FL))
#####
SRS.tp.segmented=data.frame()
for(i in 1:nrow(dates.val.srs)){
  tmp.dat=subset(grab.wq2.month.TP,Region=="SRS"&mon.cy.date==dates.val.srs$mon.cy.date[i])
  #plot(mean.TP~EuDist.frac,tmp.dat,log="y")
  mod=lm(log.TP~EuDist.frac,tmp.dat)
  change.pt=davies.test(mod,seg.Z=~EuDist.frac)
  
  x.pt=change.pt$statistic
  stat.tab=data.frame(change.pt$process)
  x.pt.stat=stat.tab[stat.tab$psi.values==x.pt,2]
  pval=change.pt$p.value
 
  tmp.rslt=data.frame(mon.cy.date=dates.val.srs$mon.cy.date[i],Region="SRS",N.sites=nrow(tmp.dat),ChangePoint=x.pt,Statistic=x.pt.stat,p.val=pval)
  SRS.tp.segmented=rbind(SRS.tp.segmented,tmp.rslt)
  print(i)
}
SRS.tp.segmented$ChangePoint.plot=with(SRS.tp.segmented,ifelse(p.val<0.05,ChangePoint,NA))

plot(ChangePoint~mon.cy.date,SRS.tp.segmented,type="l")
with(SRS.tp.segmented,points(mon.cy.date,ChangePoint.plot,pch=21))
abline(v=c(irma.landfall.FL,wilma.landfall.FL),lty=2)


test=subset(grab.wq2.month.TP,Region=="SRS"&mon.cy.date==dates.val.srs$mon.cy.date[i])
test.lm=lm(log.TP~EuDist.frac,test)

seg.lm=segmented(mblm.mod,seg.Z=~EuDist.frac)
summary(seg.lm)
plot(log.TP~EuDist.frac,test)
plot(seg.lm,add=T)


##WY
grab.wq2.WY=ddply(subset(grab.wq2,!(Station.ID%in%exclude.sites)),c("EuDist.frac","Station.ID","Region","datasource","WY"),summarise,mean.TP=mean(TP*1000,na.rm=T),N.TP=N(TP),log.TP=log(mean.TP))
grab.wq2.WY.TP=subset(grab.wq2.WY,N.TP>4)

WYS=seq(2005,2018,1)
SRS.tp.gradient.WY=data.frame()
for(i in 1:length(WYS)){
  tmp.dat=subset(grab.wq2.WY.TP,Region=="SRS"&WY==WYS[i])
  #plot(mean.TP~EuDist.frac,tmp.dat,log="y")
  mblm.mod=mblm(log.TP~EuDist.frac,tmp.dat)
  
  tmp.rslt=data.frame(WY=WYS[i],Region="SRS",N.sites=nrow(tmp.dat),slope.val=as.numeric(coefficients(mblm.mod)[2]),r2=summary.lm(mblm.mod)$r.squared,RSE=summary(mblm.mod)$sigma)
  SRS.tp.gradient.WY=rbind(SRS.tp.gradient.WY,tmp.rslt)
  print(i)
}
plot(slope.val~WY,SRS.tp.gradient.WY,type="l")
abline(v=c(WY(irma.landfall.FL),WY(wilma.landfall.FL)))



test=subset(grab.wq2.month,mon.cy.date==dates.val[4]&Region=="SRS")
plot(mean.TP~EuDist.frac,test,log="y");abline(v=0.45)
mblm(log.TP~EuDist.frac,subset(test,EuDist.frac<0.45))
mblm(log.TP~EuDist.frac,subset(test,EuDist.frac>0.45))

test.lm=lm(log(mean.TP)~EuDist.frac,test)
library(segmented)
seg.lm=segmented(test.lm,seg.Z=~EuDist.frac)
summary(seg.lm)
plot(mean.TP~EuDist.frac,test)
plot(seg.lm,add=T)
