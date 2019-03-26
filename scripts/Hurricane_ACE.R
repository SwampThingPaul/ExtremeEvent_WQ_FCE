# Weather Data2 -----------------------------------------------------------
# Data from https://www.ndbc.noaa.gov/
# calculate https://en.wikipedia.org/wiki/Accumulated_cyclone_energy
Yrs=seq(2005,2018,1)
header.val=c("YYYY", "MM", "DD", "hh", "mm", "WD", "WSPD", "GST", "WVHT","DPD", "APD", "MWD", "BAR", "ATMP", "WTMP", "DEWP", "VIS", "TIDE")

site.val=c("smkf1","nfbf1","vcaf1")
wx.dat2=data.frame()
for(j in 1:length(site.val)){
for(i in 1:length(Yrs)){
  error.val=tryCatch(read.table(paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=",site.val[j],"h",Yrs[i],".txt.gz&dir=data/historical/stdmet/"),sep="",header=T,col.names = header.val,na.strings=c("99","999","9999.0")),error=function(e) T)
  if(is.null(nrow(error.val))==T){next}else{
  tmp.dat=error.val
  tmp.dat$site=site.val[j]
  wx.dat2=rbind(tmp.dat,wx.dat2)
  print(i)}
}
print(j)
}

tmp.dat=read.table(paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=smkf1h2004.txt.gz&dir=data/historical/stdmet/"),sep="",header=T,na.strings=c("99","999","9999.0"))
names(tmp.dat)
tmp.dat$mm=as.integer(0)
tmp.dat$site="smkf1"
tmp.dat=tmp.dat[,c(header.val,"site")]
wx.dat2=rbind(wx.dat2,tmp.dat)

tmp.dat=read.table(paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=nfbf1h2004.txt.gz&dir=data/historical/stdmet/"),sep="",header=T,na.strings=c("99","999","9999.0"))
names(tmp.dat)
tmp.dat$mm=as.integer(0)
tmp.dat$site="nfbf1"
tmp.dat=tmp.dat[,c(header.val,"site")]
wx.dat2=rbind(wx.dat2,tmp.dat)

## ---

wx.dat2$DateTime.EST=with(wx.dat2,date.fun(paste(YYYY,"-",MM,"-",DD," ",hh,":",ifelse(mm==0,"00",mm),":","00",sep=""),form="%Y-%m-%d %H:%M:%S"))
wx.dat2$Date.EST=date.fun(wx.dat2$DateTime.EST)
wx.dat2$WY=WY(wx.dat2$Date.EST)

wx.dat2$WSPD[wx.dat2$WSPD>=99]=NA
wx.dat2$GST[wx.dat2$GST>=99]=NA
wx.dat2$GST=with(wx.dat2,ifelse(is.na(GST)&WSPD==0,0,GST))
wx.dat2$WSPD.knots=wx.dat2$WSPD*1.94384
wx.dat2$GST.knots=wx.dat2$GST*1.94384

wx.dat2.hr=ddply(wx.dat2,c("YYYY", "MM", "DD", "hh","Date.EST","WY"),summarise,max.GST.knt=max(GST.knots,na.rm=T))
wx.dat2.hr$hur.season=with(wx.dat2.hr,ifelse(MM%in%seq(6,11,1),"Y","N"))
wx.dat2.hr$max.GST.knt=with(wx.dat2.hr,ifelse(is.infinite(max.GST.knt)==T,0,max.GST.knt))

subset(wx.dat2,WY==2007&is.na(GST.knots)==T)
plot(subset(wx.dat2.hr,WY==2006)$max.GST.knt)

#FL_ACE=ddply(subset(wx.dat2.hr,hur.season=="Y"),"WY",summarise,N.val=N(max.GST.knt),ACE=1e-4*sum(max.GST.knt^2,na.rm=T));test
FL_ACE=ddply(subset(wx.dat2.hr,WY%in%seq(2004,2018,1)&hur.season=="Y"),"WY",summarise,N.val=N(max.GST.knt),ACE=(1e-4)*sum((max.GST.knt)^2,na.rm=T));FL_ACE
plot(ACE~WY,FL_ACE,ylim=c(60,200),ylab="ACE (x10\u2074 knots\u00B2 Yr\u207B\u00B9)",xlab="Year" )
abline(v=WY(c(irma.landfall.FL,wilma.landfall.FL)))
abline(h=c(66,111),lty=3)#below and above normal season (respectively) based on 49 year period. 
with(FL_ACE,cor.test(WY,ACE,method="spearman"))


## 
library(HURDAT)

hur.dat=get_hurdat(basin = c("AL"))
range(hur.dat$DateTime)
hur.dat$date=date.fun(hur.dat$DateTime,tz="UTC")
hur.dat$Year=as.numeric(format(hur.dat$date,"%Y"))
hur.dat$WY=WY(hur.dat$DateTime)
hur.year=ddply(hur.dat,c("Key","Name"),summarise,Year.start=min(Year,na.rm=T),WY.start=min(WY,na.rm=T),max.wind.ms=max(Wind*0.44704,na.rm=T),min.pres=min(Pressure,na.rm=T))

subset(hur.year,WY.start%in%Yrs)
