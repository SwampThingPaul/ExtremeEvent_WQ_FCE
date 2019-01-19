## Frontiers in Science 
## (Effects of Extreme Weather Events on Coastal Carbon and Nutrient Cycling)
## Version 1
## Example of Stage Elevation pre-post storm.
##
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(plyr)
library(reshape)
library(zoo);

#Custom Functions
source("D:/CommonlyUsedFunctions.r")

#the commonly used functions is also here
#source("https://github.com/SwampThingPaul/analyst_helper/blob/3bc5585bcbe2af65fdf36172dd4047e4791f594d/CommonlyUsedFunctions.r")

#Paths
setwd("D:/UF/LTER/Projects/ExtremeEvent_WQ_FCE")

paths=paste0(getwd(),c("/Exports/","/Plots/","/Data/","/GIS/"))
#Folder.Maker(paths);#One and done. Creates folders in working directory.
#Folder.Maker(paste0(getwd(),"/scripts/"))
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]



# Import Data -------------------------------------------------------------

# HydroData
dates=as.Date(c("2002-10-01","2018-10-01"))
stage.dat.P62=SFWMD.DBHYDRO.Data.daily(dates[1],dates[2],"G6150")
stage.dat.P62$Date=date.fun(stage.dat.P62$Date)
stage.dat.P62$stage.m=ft.to.m(stage.dat.P62$Data.Value)
stage.dat.P62$stage.m=ifelse(is.na(stage.dat.P62$stage.m)==T,-1,stage.dat.P62$stage.m)

wilma.period=date.fun(c("2003-10-01","2005-01-01"))
wilma.landfall.FL=date.fun("2004-10-24")

irma.period=date.fun(c("2016-10-01","2018-01-01"))
irma.landfall.FL=date.fun("2017-09-10")

ylim.val=c(0,1.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
tiff(filename=paste0(plot.path,"FigX_P62_WL.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
#png(filename=paste0(plot.path,"png/FigX_P62_WL.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1.0,0.5),oma=c(1.5,2,0.25,1));
layout(matrix(c(1:2),2,1,byrow=T),heights=c(rep(0.8,3),0.25))

xlim.val=wilma.period;xmaj=seq(xlim.val[1],xlim.val[2],"3 months");xmin=seq(xlim.val[1],xlim.val[2],"1 months")
plot(stage.m~Date,stage.dat.P62,ylim=ylim.val,xlim=xlim.val,type="n",ylab=NA,xlab=NA,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(stage.dat.P62,lines(Date,stage.m))
with(stage.dat.P62,shaded.range(Date,rep(-1,length(stage.m)),stage.m,bg="dodgerblue",lty=1))
abline(v=wilma.landfall.FL)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),0.8);axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
mtext(side=3,"Hurricane Wilma")

xlim.val=irma.period;xmaj=seq(xlim.val[1],xlim.val[2],"3 months");xmin=seq(xlim.val[1],xlim.val[2],"1 months")
plot(stage.m~Date,stage.dat.P62,ylim=ylim.val,xlim=xlim.val,type="n",ylab=NA,xlab=NA,yaxt="n",xaxt="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(stage.dat.P62,lines(Date,stage.m))
with(stage.dat.P62,shaded.range(Date,rep(-1,length(stage.m)),stage.m,bg="dodgerblue",lty=1))
abline(v=irma.landfall.FL)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"),0.8);axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
mtext(side=3,"Hurricane Irma")
mtext(side=2,"Stage Elevation (meters, NGVD29)",outer=T,line=0.5)
mtext(side=1,line=2,"Date (Month-Year)")
dev.off()