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

#GIS libraries
library(maptools)
library(classInt)
library(GISTools)
library(rgdal)
library(sp)
library(tmap)
library(tmaptools)
library(raster)
library(spatstat)
library(sf)
library(HatchedPolygons)
library(spatialEco)


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

#sample check
auto.dat2=data.frame()
sites.val=ddply(auto.dat,c("SITENAME"),summarise,N.val=N(as.numeric(DATE)))
for(i in 1:length(sites.val$SITENAME)){
  tmp=subset(auto.dat,SITENAME==sites.val$SITENAME[i])
  tmp$DATE.diff=c(NA,diff(tmp$DATE))
  auto.dat2=rbind(auto.dat2,tmp)
}

range(auto.dat2$DATE.diff,na.rm=T)
subset(auto.dat2,DATE.diff>3)

plot(TP~DATE,subset(auto.dat2,SITENAME=="SRS1d"),pch=19,col=adjustcolor("grey",0.25))
plot(TP~DATE,subset(auto.dat,SITENAME=="SRS5"))
plot(TP~DATE,subset(auto.dat,SITENAME=="SRS6"))

plot(TP~DATE,subset(auto.dat,SITENAME=="TS/PH6a"))
plot(TP~DATE,subset(auto.dat,SITENAME=="TS/PH7a"))

