
# Hurricane Index Data Analysis -------------------------------------------

# Data from https://www.ndbc.noaa.gov/
# calculate https://en.wikipedia.org/wiki/Accumulated_cyclone_energy

library(HURDAT)

hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,Year.start=min(Year,na.rm=T),WY.start=min(WY,na.rm=T),max.wind.ms=max(Wind*0.44704,na.rm=T),min.pres=min(Pressure,na.rm=T))
#subset(hur.year,WY.start%in%Yrs)
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

enp=spTransform(readOGR(GIS.path,"ENP_boundary"),utm17)
enp.buffer=gBuffer(enp,width=200*1000)
enp.buffer=SpatialPolygonsDataFrame(enp.buffer,data.frame(row.names = "buffer",width.km=200,area="ENP"))

enp.hurr=over(enp.buffer,hur.tck,returnList = T,byid=T)[[1]]
plot(enp.hurr)
enp.hurr=enp.hurr[order(enp.hurr$WY.start,enp.hurr$Key),]

hurr.N.enp=ddply(enp.hurr,c("WY.start"),summarise,N.val=N(Key))
hur.tck2=subset(hur.tck,Key%in%enp.hurr$Key)

hur.dat.enp=subset(hur.dat,Key%in%enp.hurr$Key&WY%in%seq(2005,2018,1))

plot(hur.dat.enp$Wind)

plot(Wind~DateTime,subset(hur.dat.enp,WY==2005))

max.wind=ddply(hur.dat.enp,c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
#ENP.ACE=ddply(hur.dat.enp,c("Key","Name","WY"),summarise,ACE=sum((1e-4)*sum(max(Wind,na.rm=T)^2,na.rm=T)))
ENP.ACE=ddply(max.wind,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))
ENP.ACE=rbind(ENP.ACE,data.frame(WY=c(2008,2010,2012,2015),ACE=0))

max.wind.all=ddply(subset(hur.dat,WY%in%seq(2005,2018,1)),c("WY","Key","Name"),summarise,Vmax_2=sum(Wind^2,na.rm=T))
ACE.all=ddply(max.wind.all,c("WY"),summarise,ACE=(1e-4)*sum(Vmax_2,na.rm=T))

plot(ACE~WY,ENP.ACE)

plot(ACE~WY,ACE.all,ylim=c(0,350))
with(ENP.ACE,points(WY,ACE,pch=21,bg="red"))
abline(v=WY(c(wilma.landfall.FL,irma.landfall.FL)))
