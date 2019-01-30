## Frontiers in Science 
## (Effects of Extreme Weather Events on Coastal Carbon and Nutrient Cycling)
## Version 1
## Hurricane Paths (GIS data)
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(plyr)
library(reshape)
library(zoo);
library(lubridate)

#Custom Functions
source("D:/CommonlyUsedFunctions.r")

#the commonly used functions is also here
#source("https://github.com/SwampThingPaul/analyst_helper/blob/3bc5585bcbe2af65fdf36172dd4047e4791f594d/CommonlyUsedFunctions.r")

#Paths
setwd("D:/UF/LTER/Projects/ExtremeEvent_WQ_FCE")

paths=paste0(getwd(),c("/Export/","/Plots/","/Data/","/GIS"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]

#shell.exec(system.file(getwd()))

#Helper variables 
N.mw=14.0067
P.mw=30.973762


# Water Quality Data ------------------------------------------------------
# autosampler data
na.values=c("-9999","-9999.0","-9999.00","-9999.000")
grahl.auto=read.csv(paste0(data.path,"WQ/LT_ND_Grahl_001.txt"),na.strings = na.values)
grahl.auto$DATE=date.fun(as.character(grahl.auto$DATE))
apply(grahl.auto,2,range,na.rm=T)

losada.auto=read.csv(paste0(data.path,"WQ/LT_ND_Losada_001.txt"),na.strings =na.values)
losada.auto$DATE=date.fun(as.character(losada.auto$DATE))
apply(losada.auto,2,range,na.rm=T)

rondeau.auto=read.csv(paste0(data.path,"WQ/LT_ND_Rondeau_001.txt"),na.strings =na.values)
#rondeau.auto$DATE2=date.fun(as.character(rondeau.auto$DATE))
rondeau.auto$DATE=date.fun(as.character(rondeau.auto$DATE),form="%Y-%d-%m");#Assumes YYYY-DD-MM, sent email to confirm
apply(rondeau.auto,2,range,na.rm=T)

rubio.auto=read.csv(paste0(data.path,"WQ/LT_ND_Rubio_001.txt"),na.strings =na.values)
rubio.auto$DATE=date.fun(as.character(rubio.auto$DATE))
apply(rubio.auto,2,range,na.rm=T)
rubio.auto$TP=with(rubio.auto,ifelse(TP==0.00,0.01,TP));# zero value recorded. Assume value was below detection and therefore 0.01 uM

auto.dat=rbind(grahl.auto,losada.auto,rondeau.auto,rubio.auto)
rm(grahl.auto,losada.auto,rondeau.auto,rubio.auto)
apply(auto.dat,2,range,na.rm=T)

## 
ddply(auto.dat,c("SITENAME"),summarise,min.date=min(DATE),max.date=max(DATE))

#sites=as.character(unique(auto.dat$SITENAME))
#date.fill=seq(date.fun("1996-04-07"),date.fun("2018-01-29"),"3 days")
#fill.dat=data.frame(SITENAME=sort(rep(sites,length(date.fill))),DATE=rep(date.fill,length(sites)),fill=1)
#auto.dat=merge(auto.dat,fill.dat,c("SITENAME","DATE"),all.y=T)

#ddply(auto.dat,c("SITENAME"),summarise,min.date=min(DATE),max.date=max(DATE))

alias.site=data.frame(SITENAME=c("SRS1a", "SRS1c", "SRS1d","SRS2", "SRS3", "SRS4", "SRS5", "SRS6", "TS/PH6a", "TS/PH7a","TS/PH8", "TS/Ph4", "TS/Ph5", "TS/Ph1a", "TS/Ph1b", "TS/Ph2","TS/Ph3"),
                     ALIAS=c("SRS1", "SRS1", "SRS1","SRS2", "SRS3", "SRS4", "SRS5", "SRS6", "TS/Ph6a", "TS/Ph7a","TS/Ph8", "TS/Ph4", "TS/Ph5", "TS/Ph1a", "TS/Ph1b", "TS/Ph2","TS/Ph3"))
auto.dat=merge(auto.dat,alias.site,"SITENAME")

auto.dat$TP.ugL=auto.dat$TP*P.mw
auto.dat$TN.mgL=(auto.dat$TN*N.mw)*0.001
###
pre.post.dur=2

wilma.landfall.FL=date.fun("2004-10-24")
wilma.period=c(wilma.landfall.FL-duration(pre.post.dur,"months"),wilma.landfall.FL+duration(pre.post.dur,"months"))#date.fun(c("2003-10-01","2005-01-01"))

irma.landfall.FL=date.fun("2017-09-10")
irma.period=c(irma.landfall.FL-duration(pre.post.dur,"months"),irma.landfall.FL+duration(pre.post.dur,"months"))#date.fun(c("2016-10-01","2018-01-01"))

sites=c(paste0("SRS",1:6),paste0("TS/Ph",c("1a",2,3,"6a","7a")))




# Plots
ylim.val=c(1,500);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
pch.val=rep(c(21,22,23),2)
lty.val=rep(c(1,2,3),2)
col.rmp=colorRampPalette(c("olivedrab1", "dodgerblue1"))
col=col.rmp(6)
cex.pt=1


#tiff(filename=paste0(plot.path,"/SRS_autoTP.tiff"),width=6,height=5,units="in",res=220,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,2,1.5,0.5),oma=c(3,1.5,0.5,0.25),mgp=c(3,1,0));
#layout(matrix(c(1:3),3,1,byrow=T),heights=c(1,1,0.5))
layout(matrix(c(1,2,3,3),2,2,byrow=F),widths=c(0.8,0.2))

xlim.val=wilma.period;xmaj=seq(xlim.val[1],xlim.val[2],"1 months");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
plot(TP.ugL~DATE,auto.dat,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",log="y",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:6){
  with(subset(auto.dat,ALIAS==sites[i]),pt_line(DATE,TP.ugL,lty.val[i],adjustcolor(col[i],0.5),2,pch.val[i],col[i],pt.lwd=0.1,cex=cex.pt))
}
abline(v=wilma.landfall.FL,lwd=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),1)
axis_fun(2,ymaj,ymin,ymaj,1)
mtext(side=3,"Wilma")

xlim.val=irma.period;xmaj=seq(xlim.val[1],xlim.val[2],"1 months");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
plot(TP.ugL~DATE,auto.dat,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",log="y",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:6){
  with(subset(auto.dat,ALIAS==sites[i]),pt_line(DATE,TP.ugL,lty.val[i],adjustcolor(col[i],0.5),2,pch.val[i],col[i],pt.lwd=0.1,cex=cex.pt))
}
abline(v=irma.landfall.FL)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),1)
axis_fun(2,ymaj,ymin,ymaj,1)
mtext(side=3,"Irma")

mtext(side=1,line=2,"Date (Month-Year)")
mtext(side=2,outer=T,"TP Concentration (\u03BCg L\u207B\u00B9)")
plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
leg.txt=paste0("SRS", 1:6)
legend(0.5,0.6,legend=leg.txt,pch=pch.val,pt.bg=col,col=adjustcolor(col,0.75),lty=lty.val,lwd=1.25,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(0.5,0.6,legend=leg.txt,pch=pch.val,pt.bg=col,col=c("black"),lty=0,lwd=0.5,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

ylim.val=c(0.1,5);ymaj=log.scale.fun(ylim.val[1],ylim.val[2],"major");ymin=log.scale.fun(ylim.val[1],ylim.val[2],"minor")
pch.val=rep(c(21,22,23),2)
lty.val=rep(c(1,2,3),2)
col.rmp=colorRampPalette(c("olivedrab1", "dodgerblue1"))
col=col.rmp(6)

#tiff(filename=paste0(plot.path,"/SRS_autoTN.tiff"),width=6,height=5,units="in",res=220,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,2,1.5,0.5),oma=c(3,1.5,0.5,0.25),mgp=c(3,1,0));
#layout(matrix(c(1:3),3,1,byrow=T),heights=c(1,1,0.5))
layout(matrix(c(1,2,3,3),2,2,byrow=F),widths=c(0.8,0.2))

xlim.val=wilma.period;xmaj=seq(xlim.val[1],xlim.val[2],"1 months");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
plot(TN.mgL~DATE,auto.dat,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",log="y",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:6){
  with(subset(auto.dat,ALIAS==sites[i]),pt_line(DATE,TN.mgL,lty.val[i],adjustcolor(col[i],0.5),2,pch.val[i],col[i],pt.lwd=0.1,cex=cex.pt))
}
abline(v=wilma.landfall.FL,lwd=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),1)
axis_fun(2,ymaj,ymin,ymaj,1)
mtext(side=3,"Wilma")

xlim.val=irma.period;xmaj=seq(xlim.val[1],xlim.val[2],"1 months");xmin=seq(xlim.val[1],xlim.val[2],"1 days")
plot(TN.mgL~DATE,auto.dat,type="n",ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",log="y",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:6){
  with(subset(auto.dat,ALIAS==sites[i]),pt_line(DATE,TN.mgL,lty.val[i],adjustcolor(col[i],0.5),2,pch.val[i],col[i],pt.lwd=0.1,cex=cex.pt))
}
abline(v=irma.landfall.FL)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),1)
axis_fun(2,ymaj,ymin,ymaj,1)
mtext(side=3,"Irma")

mtext(side=1,line=2,"Date (Month-Year)")
mtext(side=2,outer=T,"TN Concentration (mg L\u207B\u00B9)")
plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
leg.txt=paste0("SRS", 1:6)
legend(0.5,0.6,legend=leg.txt,pch=pch.val,pt.bg=col,col=adjustcolor(col,0.75),lty=lty.val,lwd=1.25,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(0.5,0.6,legend=leg.txt,pch=pch.val,pt.bg=col,col=c("black"),lty=0,lwd=0.5,pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()