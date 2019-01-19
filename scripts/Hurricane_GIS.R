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
GIS.path=paths[4]

GIS.path.gen="D:/_GISData"
#
nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Hurricane_TropcialStorm -------------------------------------------------
library(HURDAT)

hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,Year.start=min(Year,na.rm=T),WY.start=min(WY,na.rm=T),max.wind=max(Wind,na.rm=T),min.pres=min(Pressure,na.rm=T))

hur.dat.sp.pt=SpatialPointsDataFrame(coords=hur.dat[,c("Lon","Lat")],data=hur.dat,proj4string = CRS("+init=epsg:4269"))
#plot(hur.dat.sp.pt,pch=19)

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
sp_lines=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
#sp_lines=SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$Key))),CRS("+init=epsg:4269"))
pb=txtProgressBar(1,max=length(path),style=3)
for(i in 2:length(path)){
  tmp=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269")),data.frame(row.names=hur.id$Key,Key=hur.id$Key,Name=hur.id$Name))
  #tmp=SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269"))
  sp_lines=rbind(sp_lines,tmp)
  setTxtProgressBar(pb,i)
}
chk=data.frame(gIsValid(sp_lines,byid=T))
chk$Key=rownames(chk)
colnames(chk)=c("geo.valid","Key")
subset(chk,geo.valid=="FALSE")

hur.tck=spTransform(sp_lines,utm17)
hur.tck=merge(hur.tck,hur.year,by.x=c("Key","Name"),by.y=c("Key","Name"))

subset(hur.tck,Key%in%c("AL112017","AL252005"))
plot(subset(hur.tck,Key%in%c("AL112017","AL252005")))
#writeOGR(subset(hur.tck,Key%in%c("AL112017","AL252005")),GIS.path,"Wilma_Irma_bestrack",driver="ESRI Shapefile")
Irma_Wilma=subset(hur.tck,Key%in%c("AL112017","AL252005"))
Irma_Wilma$Name=factor(Irma_Wilma$Name,levels=as.character(unique(Irma_Wilma$Name)))


date.fun("2005-10-24 10:30:00",form="%F %X",tz="UTC")

land.fall=data.frame(Storm=c("Wilma","Irma"),
                     Notes=c("Florida landfall - near Cape Ramono","2nd Florida landfall - near Marco island"),
                     Date_Time_UTC=c(date_fun("2005-10-24 10:30:00",form="%F %X",tz="UTC"),date_fun("2017-09-10 19:30:00",form="%F %X",tz="UTC")),
                     Lat=c(25.9,25.9),
                     Lon=c(81.7,81.7),
                     Pressure.mb=c(950,936),
                     WindSp.kt=c(105,100))
subset(sp_lines,Name%in%c("IRMA","WILMA"))@data


shoreline=readOGR(GIS.path.gen,"FloridaShoreline")
shoreline=spTransform(shoreline,utm17)

plot(subset(hur.tck,Key%in%c("AL112017","AL252005")))
plot(shoreline)
plot(subset(hur.tck,Key%in%c("AL112017","AL252005")),add=T)

tmap_mode("view")
tm_shape(Irma_Wilma)+tm_lines(col="Name",lwd=5)
