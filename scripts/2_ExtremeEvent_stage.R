##

# Water Level -------------------------------------------------------------
#stg.sites=c("ENPSR","ENPBR","ENPCN","NP-P38","NP-P34","NP-P36","NP-A13","NP-205","NP-202","G-3273","NP-201","NESRS2")
stg.sites.SRS=c("ENPSR","ENPGI","ENPCN","NP-P36","NP-202","NP-201","NESRS2");#SRS
stg.sites.TS=c("NP-TSB","R-127","NP-TSH","NP-146","ENPLM")#,"ENPTR");#TS
stg.sites=data.frame(STATION=c(stg.sites.SRS,stg.sites.TS),Region=c(rep("SRS",length(stg.sites.SRS)),rep("TS",length(stg.sites.TS))))
sfwmd.mon.stage.trans=subset(sfwmd.mon,ACTIVITY_S%in%c("Stage","Well")&STATUS=="Active"&STATION%in%c(as.character(stg.sites$STATION),"NP-31W","S12A_T","S12B_T","S12C_T","S12D_T","S333_T","S334_H"))
sfwmd.mon.q=subset(sfwmd.mon,ACTIVITY_S%in%c("Flow")&STATUS=="Active")
#paste(sfwmd.mon.stage.trans@data$STATION,collapse = "/")
#tm_shape(sfwmd.mon.stage.trans)+tm_dots()

## Euclidean Distance from GOM/FB
GOMr=raster(xmn=extent(sfwmd.mon)[1],xmx=extent(sfwmd.mon)[2],ymn=extent(sfwmd.mon)[3],ymx=extent(sfwmd.mon)[4],crs=utm17)
GOMr=setValues(GOMr,0)
GOMr=mask(GOMr,subset(sfwmd.mon.stage.trans,STATION=="ENPSR"))
GOMrD=distance(GOMr)
#plot(GOMrD);plot(sfwmd.mon.stage.trans,add=T,pch=19)
stg.srs=subset(sfwmd.mon.stage.trans,ACTIVITY_S%in%c("Stage","Well")&STATUS=="Active"&STATION%in%c(as.character(subset(stg.sites,Region=="SRS")$STATION),"S12A_T","S12B_T","S12C_T","S12D_T","S333_T","S334_H"))
stg.sites.eudist.srs=data.frame(STATION=stg.srs$STATION,EuDist.m=raster::extract(GOMrD,stg.srs))
stg.sites.eudist.srs$Frac.dist=1-stg.sites.eudist.srs$EuDist.m/max(stg.sites.eudist.srs$EuDist.m)

FBr=raster(xmn=extent(sfwmd.mon)[1],xmx=extent(sfwmd.mon)[2],ymn=extent(sfwmd.mon)[3],ymx=extent(sfwmd.mon)[4],crs=utm17)
FBr=setValues(FBr,0)
FBr=mask(FBr,subset(sfwmd.mon.stage.trans,STATION=="ENPLM"))
FBrD=distance(FBr)
#plot(FBrD);plot(sfwmd.mon.stage.trans,add=T,pch=19)
stg.ts=subset(sfwmd.mon.stage.trans,ACTIVITY_S%in%c("Stage","Well")&STATUS=="Active"&STATION%in%c(as.character(subset(stg.sites,Region=="TS")$STATION),"NP-31W"))
stg.sites.eudist.ts=data.frame(STATION=stg.ts$STATION,EuDist.m=raster::extract(FBrD,stg.ts))
stg.sites.eudist.ts$Frac.dist=1-stg.sites.eudist.ts$EuDist.m/max(stg.sites.eudist.ts$EuDist.m)
stg.sites=merge(stg.sites,rbind(stg.sites.eudist.srs,stg.sites.eudist.ts),"STATION")

#WL.ts=openxlsx::read.xlsx(paste0(data.path,"/WL_DBKeys_REPORT.xlsx"))
#WL.ts$Start.Date=date.fun(openxlsx::convertToDate(WL.ts$Start.Date))
#WL.ts$End.Date=date.fun(openxlsx::convertToDate(WL.ts$End.Date))

library(ncdf4)
dem=nc_open(paste0(data.path,"EDEN/eden_dem_cm_oc11.nc"))
lon <- ncvar_get(dem,"x")
nlon <- dim(lon)
lat <- ncvar_get(dem,"y")
nlat <- dim(lat)
dname="dem"
tmp_array <- ncvar_get(dem,dname)
dlname <- ncatt_get(dem,dname,"long_name")
dunits <- ncatt_get(dem,dname,"units")
fillvalue <- ncatt_get(dem,dname,"_FillValue")
title <- ncatt_get(dem,0,"title")
institution <- ncatt_get(dem,0,"institution")
datasource <- ncatt_get(dem,0,"source")
references <- ncatt_get(dem,0,"references")
history <- ncatt_get(dem,0,"history")
Conventions <- ncatt_get(dem,0,"Conventions")
tmp_array[tmp_array==fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,1])))
tmp_slice <- tmp_array
# Map as raster
dem.r <- raster(t(tmp_slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=utm +zone=17 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
dem.r <- flip(dem.r, direction='y')

dem.r.FCE=mask(dem.r,enp.shoreline)
dev.off();plot(dem.r.FCE)

select.stg=subset(sfwmd.mon,ACTIVITY_S%in%c("Stage","Well")&STATUS=="Active"&STATION%in%stg.sites$STATION)
tm_shape(enp.shoreline)+tm_borders(col="grey10",lwd=1,alpha=0.5)+
  tm_shape(enp.slough.clip)+tm_polygons("grey",alpha=0.5)+
  tm_shape(fce.sites)+tm_dots(col="grey",size=0.2)+
  tm_shape(select.stg)+tm_dots(col="dodgerblue1",size=0.06)

stage.sites=data.frame(Station=c("ENPSR","ENPCN","ENPGI","NESRS2","NP-146","NP-201","NP-202","NP-P36","NP-TSB","NP-TSH","ENPLM"),
                       DBKEY=c("63681","63599","63623","01218","H2428","06719","06720","06718","H2442","07090","63644"))
stage.sites=merge(stage.sites,stg.sites,by.x="Station",by.y="STATION")
stg.elev.dat=data.frame(Station=c("ENPCN","ENPGI","ENPLM","ENPSR","NESRS2","NP-146","NP-201","NP-202","NP-P36","NP-TSB","NP-TSH"),
                        GndElev.ft.ngvd29=c(1.11,-1.48,0.28,-1.92,5.85,0.70,6.90,5.86,3.50,3.83,1.43))
stage.sites=merge(stage.sites,stg.elev.dat,"Station")
#paste(stage.sites$DBKEY,collapse="/")

#tmap_mode("view")
bbox.val=bbox(enp.slough.clip)
stage.map=tm_shape(enp.shoreline)+tm_borders(col="grey")+
  tm_shape(dem.r.FCE,bbox=bbox(enp.shoreline))+tm_raster(palette="-Greys",n=6,alpha=0.5,title="Elevation\n(cm, NGVD88)")+
  tm_shape(enp.slough.clip)+tm_borders(alpha=0.5)+
  tm_shape(select.stg)+tm_symbols(col="white",size=0.75)+tm_text("STATION",just="bottom",col="red",size=0.5)+
  tm_add_legend(type="symbol",label="Stage Monitoring",col="white",shape=21)+
  tm_layout(legend.show = T, legend.outside = T)
#stage.map
#tmap_save(stage.map,paste0(plot.path,"stage_monitoring.png"),width = 6.5,height=4,units="in",dpi=200)


# Stage Data --------------------------------------------------------------
sdate=as.Date("2004-05-01");# Begining of Florida WY2005
edate=as.Date("2018-04-30");# End of Florida WY2018

#paste(stage.sites$DBKEY,collapse="/")
WL.dat=data.frame()
for(i in 1:nrow(stage.sites)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,stage.sites$DBKEY[i])
  tmp$DBKEY=as.character(stage.sites$DBKEY[i])
  WL.dat=rbind(WL.dat,tmp)
  print(i)
}
WL.dat=merge(WL.dat,stage.sites,c("DBKEY","Station"))
WL.dat$Date=date.fun(WL.dat$Date)
WL.dat$Data.Value=with(WL.dat,ifelse(Data.Value==0,NA,Data.Value))
WL.dat$Data.Value.m=ft.to.m(WL.dat$Data.Value)
WL.dat$depth.cm=with(WL.dat,ft.to.m(Data.Value-GndElev.ft.ngvd29)*100)
WL.dat$WL_LTgnd=with(WL.dat,ifelse(depth.cm<0,0,1))
WL.dat$hydroperiod.90dfreq=with(WL.dat,ave(WL_LTgnd,Station,FUN=function(x)c(rep(NA,89),rollsum(x,90)/90)))
plot(hydroperiod.90dfreq~Date,subset(WL.dat,Station=="NP-202"))

c("ENPCN", "ENPGI", "ENPLM", "ENPSR", "NESRS2", "NP-146", "NP-201","NP-202", "NP-P36", "NP-TSB", "NP-TSH")

WL.dat=merge(WL.dat,data.frame(Station=c("NP-202", "NP-P36","NP-P36","ENPCN","ENPGI","ENPSR","NP-TSB","NP-146","ENPLM"),
                               WQSite=c(paste0("SRS",c("1d",2:6)),paste0("TS/PH",c(2,3,"7a")))),"Station",all.x=T)

WL.dat.xtab=cast(subset(WL.dat,Station%in%c("NP-202","ENPGI","ENPSR","NP-TSB","ENPLM")),Date~Station,value="Data.Value",fun.aggregate=function(x) mean(ft.to.m(x),na.rm=T))
WL.dat.xtab=rename(WL.dat.xtab,c("NP-202"="NP202","NP-TSB"="NPTSB"))
WL.dat.xtab$SRS=with(WL.dat.xtab,(ENPSR-NP202)/subset(stage.sites,Station=="NP-202")$EuDist.m)
WL.dat.xtab$TS=with(WL.dat.xtab,(ENPLM-NPTSB)/subset(stage.sites,Station=="NP-TSB")$EuDist.m)

WL.dat.xtab.melt=melt(data.frame(WL.dat.xtab)[,c("Date","SRS","TS")],id.vars="Date",variable_name="Region")
WL.dat.xtab.melt=rename(WL.dat.xtab.melt,c("value"="gradient"))

WL.dat2=merge(WL.dat,analysis.periods,by.x="Date",by.y="DATE")
WL.dat2$PeriodStorm=with(WL.dat2,paste(Period,Storm,sep="_"))

storm.order=data.frame(PeriodStorm=c("Pre_Wilma","Post_Wilma","Pre_Irma","Post_Irma"),PeriodStorm2=c("1Pre_Wilma","2Post_Wilma","3Pre_Irma","4Post_Irma"))
WL.dat2=merge(WL.dat2,storm.order,"PeriodStorm")
WL.dat2$Period=factor(WL.dat2$Period,levels=c("Pre","Post"))
WL.dat2$Storm=factor(WL.dat2$Storm,levels=c("Wilma","Irma"))
WL.dat2$Station_Storm=with(WL.dat2,paste(Station,Storm,sep="_"))

stage.sites=stage.sites[order(stage.sites$Region,stage.sites$Frac.dist),]
stage.sites=subset(stage.sites,!(Station%in%c("NESRS2","NP-201")))

sample.size=ddply(WL.dat2,c("Station","Storm"),summarise,N.val=N(Data.Value))
sample.size=subset(sample.size,N.val>61)
sample.size$Station_Storm=with(sample.size,paste(Station,Storm,sep="_"))
ddply(subset(WL.dat2,Station_Storm%in%sample.size$Station_Storm),c("Station","Storm"),summarise,N.val=N(Data.Value))
#sample.size$PeriodStorm=with(sample.size,paste(Period,Storm,sep="_"))
kw.test.storms=ddply(subset(WL.dat2,Station_Storm%in%sample.size$Station_Storm),c("Station","Storm"),summarise,chisq=kruskal.test(Data.Value,as.factor(Period))$statistic,pval=kruskal.test(Data.Value,as.factor(Period))$p.value)
#kw.test.storms=ddply(subset(WL.dat2,Station_Storm%in%sample.size$Station_Storm),c("Station","Storm"),summarise,chisq=kruskal.test(depth.cm,as.factor(Period))$statistic,pval=kruskal.test(depth.cm,as.factor(Period))$p.value)
kw.test.storms$Chisq.pval=with(kw.test.storms,paste(format(round(chisq,1),nsmall=1),ifelse(pval<0.01,"<0.01",round(pval,2)),sep=", "))
kw.test.storms=merge(kw.test.storms,stage.sites[,c("Station","Region","Frac.dist")],"Station")

kw.test.storms=kw.test.storms[order(kw.test.storms$Region,kw.test.storms$Frac.dist),]
kw.test.storms.rslt=spread(kw.test.storms[,c("Region","Frac.dist","Station","Storm","Chisq.pval")],Storm,Chisq.pval)
#write.csv(kw.test.storms.rslt[,c("Region","Station","Wilma","Irma")],paste0(export.path,"kw_WLCompare.csv"),row.names=F)

plot(Data.Value.m~Date,subset(WL.dat2,Station=="NP-P36"&Storm=="Irma"))

region.val=c("SRS","TS")
par(family="serif",mar=c(1.5,3.2,0.5,0.5),oma=c(2,1.5,1,0.25),mgp=c(3,1,0));
layout(matrix(1:10,5,2,byrow=F))
for(j in 1:2){
  sites.val=subset(stage.sites,Region==region.val[j])$Station
  for(i in 1:length(sites.val)){
    boxplot(depth.cm~Period+Storm,subset(WL.dat2,Station==as.character(sites.val[i])))
  }
}

#tiff(filename=paste0(plot.path,"FigX_PrePostStage.tiff"),width=5,height=6.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,3.2,0.5,0.5),oma=c(2,1.5,1,0.25),mgp=c(3,1,0));
layout(matrix(1:10,5,2,byrow=F))
region.val=c("SRS","TS")
for(j in 1:2){
  sites.val=subset(stage.sites,Region==region.val[j])$Station
  ylim.min=if(j==1){c(2.3,1.4,0.6,0.3,-0.25)}else{c(1.2,0.6,0.3,0)}
  #ylim.min=if(j==1){c(7.5,4.5,2,1,-1)}else{c(4,2,1,0)}
  ylim.max=if(j==1){c(3,2,1.2,1.2,0.8)}else{c(2,1.2,1,1)}
  #ylim.max=if(j==1){c(9.75,6.5,4,4,2.5)}else{c(7,4,3,3)}
  by.y=if(j==1){c(0.25,0.25,0.25,0.25,0.25)}else{c(0.25,0.25,0.25,0.25)}
  for(i in 1:length(sites.val)){
    ylim.val=c(ylim.min[i],ylim.max[i]);ymaj=c(seq(ylim.val[1],ylim.val[2],by.y[i]));ymin=seq(ylim.val[1],ylim.val[2],by.y[i]/2)
    x=boxplot(Data.Value.m~Period+Storm,subset(WL.dat2,Station==as.character(sites.val[i])),outline=F,col=c("grey","white"),ylim=ylim.val,yaxt="n",xaxt="n",lwd=0.25)
    abline(v=2.5,lty=3,col="grey")
    #abline(h=subset(stage.sites,Region==region.val[j]&Station==as.character(sites.val[i]))$GndElev.ft.ngvd29,lty=2)
    axis_fun(2,ymaj,ymin,format(ymaj,nsmall = 1),1)
    if(j==1&i==5){axis_fun(1,line=-0.5,c(1:4),c(1:4),rep(c("Pre","Post"),2),1)}else(axis_fun(1,1:4,1:4,NA,1))
    if(j==2&i==4){axis_fun(1,line=-0.5,c(1:4),c(1:4),rep(c("Pre","Post"),2),1)}else(axis_fun(1,1:4,1:4,NA,1))
    if(j==1&i==5){axis(1,line=0.5,lty=0,c(1.5,3.5),c("Wilma","Irma"))}
    if(j==2&i==4){axis(1,line=0.5,lty=0,c(1.5,3.5),c("Wilma","Irma"))}
    if(i==1){mtext(side=3,region.val[j])}
    mtext(side=3,line=-1.1,as.character(sites.val[i]),cex=0.8)
    #DT=with(subset(WL.dat2,Station==as.character(sites.val[i])),dunn.test(Data.Value,PeriodStorm2))
    #DT2=merge(cldList(P.adjusted ~ comparison,data=DT,threshold = 0.05),data.frame(PeriodStorm2=storm.order$PeriodStorm2),by.x="Group",by.y="PeriodStorm2",all.y=T)
    #DT2=DT2[order(as.numeric(substring(DT2$Group,1,1))),]
    #DT2[is.na(DT2)]=" "
    #ltr=toupper(DT2$Letter)
    #ltr[is.na(ltr)]=" "
    #dunn.letters(4,1:4,x$stats[5,],ltr,"red",1)
    mtext(side=2,line=-0.5,outer=T,"Stage Elevation (meters, NGVD29)")
  }
}
dev.off()

CA365=SFWMD.DBHYDRO.Data.daily(sdate,edate,"16538")
CA365$Date=date.fun(CA365$Date)
CA365=merge(CA365,analysis.periods,by.x="Date",by.y="DATE")
CA365$Period=factor(CA365$Period,levels=c("Pre","Post"))
CA365$Storm=factor(CA365$Storm,levels=c("Wilma","Irma"))
CA365$Station_Storm=with(CA365,paste(Station,Storm,sep="_"))
CA365$PeriodStorm=with(CA365,paste(Period,Storm,sep="_"))
boxplot(Data.Value~Period+Storm,CA365)


WL.dat2$depth.ft=with(WL.dat2,Data.Value-GndElev.ft.ngvd29)
WL.dat2$DateTime.EST=with(WL.dat2,date.fun(paste(Date,"01:00:00"),form="%Y-%m-%d %H:%M:%S"))#fake Date Time

plot(depth.ft~Date,subset(WL.dat2,Station=="NP-201"&Storm=="Wilma"))
WL.dat2.xtab=cast(subset(WL.dat2,Station%in%c("NP-202","ENPGI","ENPSR","NP-TSB","ENPLM")),Date+Period+Storm+PeriodStorm~Station,value="Data.Value",fun.aggregate=function(x) mean(ft.to.m(x),na.rm=T))
WL.dat2.xtab=rename(WL.dat2.xtab,c("NP-202"="NP202","NP-TSB"="NPTSB"))
WL.dat2.xtab$SRS.gradient=with(WL.dat2.xtab,(ENPSR-NP202)/subset(stage.sites,Station=="NP-202")$EuDist.m)
WL.dat2.xtab$TS.gradient=with(WL.dat2.xtab,(ENPLM-NPTSB)/subset(stage.sites,Station=="NP-TSB")$EuDist.m)
WL.dat2.xtab$storm.day=with(WL.dat2.xtab,ifelse(Storm=="Wilma",difftime(Date,wilma.landfall.FL,units="days"),difftime(Date,irma.landfall.FL,units="days")))
rf.dat2$DateTime.EST=with(rf.dat2,date.fun(paste(Date.EST,"01:00:00"),form="%Y-%m-%d %H:%M:%S"))#fake Date Time
(subset(stage.sites,Station=="ENPSR")$GndElev.ft.ngvd29-subset(stage.sites,Station=="NP-202")$GndElev.ft.ngvd29)/subset(stage.sites,Station=="NP-202")$EuDist.m
(subset(stage.sites,Station=="ENPLM")$GndElev.ft.ngvd29-subset(stage.sites,Station=="NP-TSB")$GndElev.ft.ngvd29)/subset(stage.sites,Station=="NP-TSB")$EuDist.m

plot(SRS.gradient~Date,subset(WL.dat2.xtab,Storm=="Wilma"))
plot(SRS.gradient~Date,subset(WL.dat2.xtab,Storm=="Irma"))
plot(SRS.gradient~storm.day,subset(WL.dat2.xtab,Storm=="Wilma"))
with(subset(WL.dat2.xtab,Storm=="Irma"),points(storm.day,SRS.gradient,pch=21,bg="grey"))

plot(TS.gradient~Date,subset(WL.dat2.xtab,Storm=="Wilma"))
plot(TS.gradient~Date,subset(WL.dat2.xtab,Storm=="Irma"))
plot(TS.gradient~storm.day,subset(WL.dat2.xtab,Storm=="Wilma"))
with(subset(WL.dat2.xtab,Storm=="Irma"),points(storm.day,TS.gradient,pch=21,bg="grey"))

##

ddply(subset(WL.dat2.xtab,storm.day<0),"Storm",summarise,min.val=min(SRS.gradient,na.rm=T),max.val=max(SRS.gradient,na.rm=T))
ddply(subset(WL.dat2.xtab,storm.day<0),"Storm",summarise,min.val=min(TS.gradient,na.rm=T),max.val=max(TS.gradient,na.rm=T))

kruskal.test(SRS.gradient~Storm,subset(WL.dat2.xtab,storm.day<0))
boxplot(SRS.gradient~Storm,subset(WL.dat2.xtab,storm.day<0));#add to supplemental material

subset(WL.dat2.xtab,storm.day%in%c(-15:0)&Storm=="Irma")

kruskal.test(TS.gradient~Storm,subset(WL.dat2.xtab,storm.day<0))
boxplot(TS.gradient~Storm,subset(WL.dat2.xtab,storm.day<0))


xlim.val=c(-60,60);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
#tiff(filename=paste0(plot.path,"FigX_stage_gradient.tiff"),width=6.5,height=4.75,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1.5,3.2,0.25,0.25),oma=c(1,1.5,1,0.25),mgp=c(3,1,0));
layout(matrix(c(1:5,5),3,2,byrow=T),heights=c(0.2,0.8,0.2))

ylim.val=c(-0.5,3);by.y=2;alt.ylim.min=0;ymaj=seq(alt.ylim.min,ylim.val[2],by.y);ymin=seq(alt.ylim.min,ylim.val[2],by.y/2)
plot(NP202~storm.day,WL.dat2.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(WL.dat2.xtab,Storm=="Wilma"),lines(storm.day,NP202,lwd=1.25,lty=1))
with(subset(WL.dat2.xtab,Storm=="Wilma"),lines(storm.day,ENPSR,lwd=1.25,lty=2))
with(subset(WL.dat2.xtab,Storm=="Irma"),lines(storm.day,NP202,lwd=1.25,lty=1,col="indianred1"))
with(subset(WL.dat2.xtab,Storm=="Irma"),lines(storm.day,ENPSR,lwd=1.25,lty=2,col="indianred1"))
axis_fun(1,xmaj,xmin,NA,1)
axis_fun(2,ymaj,ymin,ymaj,1)
mtext(side=3,"SRS")
mtext(side=2,line=2,"Stage\n(m, NGVD29)",cex=0.75)

ylim.val=c(0,2);by.y=1;alt.ylim.min=0;ymaj=seq(alt.ylim.min,ylim.val[2],by.y);ymin=seq(alt.ylim.min,ylim.val[2],by.y/2)
plot(NPTSB~storm.day,WL.dat2.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(WL.dat2.xtab,Storm=="Wilma"),lines(storm.day,NPTSB,lwd=1.25,lty=1))
with(subset(WL.dat2.xtab,Storm=="Wilma"),lines(storm.day,ENPLM,lwd=1.25,lty=2))
with(subset(WL.dat2.xtab,Storm=="Irma"),lines(storm.day,NPTSB,lwd=1.25,lty=1,col="indianred1"))
with(subset(WL.dat2.xtab,Storm=="Irma"),lines(storm.day,ENPLM,lwd=1.25,lty=2,col="indianred1"))
axis_fun(1,xmaj,xmin,NA,1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=3,"TS")

ylim.val=c(-6.0e-5,-3.9e-5);by.y=0.5e-5;alt.ylim.min=-4e-5;ymaj=seq(ylim.val[1],alt.ylim.min,by.y);ymin=seq(ylim.val[1],alt.ylim.min,by.y/2)
plot(SRS.gradient~storm.day,WL.dat2.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(subset(WL.dat2.xtab,Storm=="Wilma"),points(storm.day,SRS.gradient,pch=21,bg="grey"))
with(subset(WL.dat2.xtab,Storm=="Wilma"),pt_line(storm.day,SRS.gradient,2,"grey",1,21,"grey"))
#with(subset(WL.dat2.xtab,Storm=="Irma"),points(storm.day,SRS.gradient,pch=19))
with(subset(WL.dat2.xtab,Storm=="Irma"),pt_line(storm.day,SRS.gradient,2,"black",1,19,"black"))
axis_fun(1,xmaj,xmin,xmaj,1)
axis_fun(2,ymaj,ymin,format(ymaj/1e-5),1)
mtext(side=2,line=2.25,"Surface Water Gradient (x10\u207B\u2075; m m\u207B\u00B9)")
#mtext(side=1,line=2,"Days till storm landfall")

ylim.val=c(-5.5e-5,-2e-5);by.y=1e-5;alt.ylim.min=-2.5e-5;ymaj=seq(ylim.val[1],alt.ylim.min,by.y);ymin=seq(ylim.val[1],alt.ylim.min,by.y/2)
plot(TS.gradient~storm.day,WL.dat2.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(subset(WL.dat2.xtab,Storm=="Wilma"),points(storm.day,SRS.gradient,pch=21,bg="grey"))
with(subset(WL.dat2.xtab,Storm=="Wilma"),pt_line(storm.day,TS.gradient,2,"grey",1,21,"grey"))
#with(subset(WL.dat2.xtab,Storm=="Irma"),points(storm.day,SRS.gradient,pch=19))
with(subset(WL.dat2.xtab,Storm=="Irma"),pt_line(storm.day,TS.gradient,2,"black",1,19,"black"))
axis_fun(1,xmaj,xmin,xmaj,1)
axis_fun(2,ymaj,ymin,format(ymaj/1e-5),1)

plot(0:1,0:1,ylab=NA,xlab=NA,axes=F,type="n")
mtext(side=3,line=-1.5,"Days till storm landfall")
legend.text=c("Upstream Stage (Wilma)","Downstream Stage (Wilma)","Upstream Stage (Irma)","Downstream Stage (Irma)","Wilma","Irma")
ln.col=c(rep("black",2),rep("indianred1",2),"grey","black")
pt.cols=c(NA,NA,NA,NA,"grey","black")
legend(0.5,0.5,legend=legend.text,lty=c(1,2,1,2,2,2),col=c(rep("black",2),rep("indianred1",2),"grey","black"),lwd=c(1,1,1,1,1),pch =c(NA,NA,NA,NA,21,19),pt.bg=pt.cols,pt.cex=1,ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(0.5,0.5,legend=legend.text,lty=c(1,2,1,2,NA,NA),col=c(rep("black",2),rep("indianred1",2),"black","black"),lwd=c(1,1,1,1,1),pch =c(NA,NA,NA,NA,21,19),pt.bg=pt.cols,pt.cex=1,ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

