
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
library(HURDAT)

library(dunn.test)
library(rcompanion)
library(smatr)
library(hydrostats)
library(tidyr)
library(gvlma)
library(mblm)
#library(trend)
#library(segmented)
#library(vegan)

#GIS libraries
library(tmap)
library(tmaptools)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)

#For disturbance CSF
library(MuMIn)
library(car)
library(nlme)
library(lattice)
library(visreg)
library(piecewiseSEM)

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

enp=spTransform(readOGR(GIS.path,"ENP_boundary"),utm17)
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

# Water Quality Analysis #1 -----------------------------------------------
#source(paste0(getwd(),"/scripts/5_ExtremeEvent_WQEval.r"),echo=T)