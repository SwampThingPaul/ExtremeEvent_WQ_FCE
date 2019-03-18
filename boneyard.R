
## WQ data viz 
library(OceanView)
Cross2 <- db2cross(data.frame(subset(grab.wq2.month,Region=="SRS")), row = "dec.year", 
                   col = "EuDist.frac", val = "mean.TP")

image2D(x = Cross2$x, y = Cross2$y, z = Cross2$z, 
        ylim = c(0, 1.01), resfac=1,log = "", 
        xlab = "Date", ylab = "Distance Downstream", clab = "TP \u03BCg L\u207B\u00B9",NAcol="grey20",yaxt="n")
axis_fun(2,line=-1.75,axisLine=-1.75,seq(0,1,0.2),seq(0,1,0.1),format(seq(0,1,0.2)))
abline(v=decimal_date(c(irma.landfall.FL,wilma.landfall.FL)),col="white",lty=2)

library(akima)
im=with(subset(grab.wq2.month,Region=="SRS"),interp(dec.year,EuDist.frac,mean.TP,duplicate="mean"))
with(im,image(x,y,z,col=grey.colors(10)))

image2D(im$x,im$y,im$z)



###
grab.wq2$log.TP=log(grab.wq2$TP)
plot(mean.TP~EuDist.frac,subset(grab.wq2.month,Region=="SRS"),log="y")

fce.wq=merge(fce.wq,wq.sites.sp@data[,c("STATION","datasource","EuDist.frac")],by.x="Station.ID",by.y="STATION",all.x=T)

date.srs=unique(sort(subset(grab.wq2.month,Region=="SRS")$mon.cy.date))
uptake.val=data.frame()
for(i in 1:length(date.srs)){
  mod=lm(log(mean.TP)~EuDist.frac,subset(grab.wq2.month,Region=="SRS"&mon.cy.date==date.srs[i]))
  kval.llm=as.numeric(coefficients(mod)[2])*-1
  uptake.val=rbind(uptake.val,data.frame(mon.cy.date=date.srs[i],kval.llm=kval.llm))
  print(paste(i,"/",length(date.srs)))
}
uptake.val$Sw=1/uptake.val$kval.llm
plot(Sw~mon.cy.date,uptake.val)
plot(mean.TP~EuDist.frac,subset(grab.wq2.month,datasource=="FCE"&Region=="SRS"&mon.cy.date==date.srs[2]&is.na(mean.TP)==F))
abline(v=c(irma.landfall.FL,wilma.landfall.FL))




kCstar0.fun=function(pars,xdata,ydata){
  C0 = pars[1]
  kval = pars[2]
  ypred = C0*exp(kval*xdata)
  loglike = -0.5*length(ydata)*log(2*pi)-sum((ypred-ydata)^2)#maximize likelihood (minimize negative)
  return(-loglike)
}
uptake.val.optim=data.frame()
for(i in 1:length(date.srs)){
  tmp.dat=subset(grab.wq2.month,Region=="SRS"&mon.cy.date==date.srs[i]&is.na(mean.TP)==F)
  y=tmp.dat$mean.TP
  x=tmp.dat$EuDist.frac
  plot(y~x)
  params=c(10,2)
  optout=optim(par=params,fn=kCstar0.fun,xdata=x,ydata=y,method="L-BFGS-B",lower=c(2,2,0),upper=c(100,1000,20))
  opt.mod=optout$par[1]+(optout$par[2]-optout$par[1])*exp(-optout$par[3]*x)
  y.mod=optout$par[1]*exp(-optout$par[2]*x)
  RSS=sum((y-y.mod)^2)# Residual Sum of Squares
  TSS=sum((y-mean(y))^2)#Total Sum of Squares
  R2=1-(RSS/TSS) #R2
  RMSE=sqrt(mean((y.mod-y)^2))# Root mean standard error
  tmp=data.frame(mon.cy.date=date.srs[i],N.vals=length(y),loglike=-optout$value,Cstar =0, C0=optout$par[1],kval=optout$par[2],RSS=RSS,TSS=TSS,R2=R2,RMSE=RMSE)
  uptake.val.optim=rbind(tmp,uptake.val.optim)
  print(paste(i,"/",length(date.srs)))
}
uptake.val.optim


# Load evaluation  --------------------------------------------------------

test=lm(TOCLoad~WY,subset(coastal.wy.region,Region=="SRS"))
library(segmented)
test.seg=segmented(test,seg.Z=~WY)
summary(test.seg)
plot(TOCLoad~WY,subset(coastal.wy.region,Region=="SRS"))
plot(test.seg,add=T)

layout(matrix(1:4,2,2));plot(test)
dev.off()
plot(test,which=1)
plot(test,which=2)
plot(test,which=3)
plot(test,which=4)
plot(test,which=5)
plot(test,which=6)

plot(test.seg$fitted.values,test.seg$residuals);#Plot 1
with(lowess(test.seg$fitted.values,test.seg$residuals),lines(x,y,col="blue"))

qqnorm(test.seg$residuals);qqline(test.seg$residuals,lty=3);#Plot 2

plot(test$fitted.values,rstudent(test.seg))








































high.spells.v2(subset(HW.flow.wq.da.region,Region=="SRS"&WY==WY(irma.landfall.FL))[,c("Date","Q")],WY.type="FL")

high.spells.v2(subset(HW.flow.wq.da.region,Region=="SRS"&Date%in%c(analysis.period.pre.wilma,analysis.period.post.wilma))[,c("Date","Q")],WY.type="FL")

plot(TFlow~Date.EST,subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(wilma.landfall.FL)))
abline(v=wilma.landfall.FL)
plot(TFlow~Date.EST,subset(HW.flow.wq.da.region,Region=="TSPh"&WY==WY(irma.landfall.FL)))
abline(v=irma.landfall.FL)

plot(TFlow~Date.EST,subset(HW.flow.wq.da.region,Region=="SRS"&WY==2018))
with(subset(HW.flow.wq.da.region,Region=="SRS"&TFlow==max(subset(HW.flow.wq.da.region,Region=="SRS")$TFlow,na.rm=T)),points(Date.EST,TFlow,pch=21,bg="red"))
abline(h=quantile(subset(HW.flow.wq.da.region,Region=="SRS")$TFlow,probs=c(0.9)),col="Red")
abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="black",lty=2)









## Coastal discharge
#tmap_mode("view")
tm_shape(enp.slough.clip)+tm_polygons(alpha=0.5)+
  tm_shape(sfwmd.mon.q)+tm_dots(size=0.05)

coastal.flow.sites=data.frame(Station=c("HARNEYFLAM","SHARKRIVBG","02290878","TAYLORS3","EASTCK","MUD_CRKM"),
                              RiverName=c("Harney","Shark","Broad","Taylor","East","Mud"),
                              Region=c(rep("SRS",3),rep("TS",3)),
                              DBKEY=c("90144","90243","89871","AN865","AM250","AN862"))

coastal.q=data.frame()
for(i in 1:nrow(coastal.flow.sites)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,coastal.flow.sites$DBKEY[i])
  tmp$DBKEY=as.character(coastal.flow.sites$DBKEY[i])
  tmp$Station=as.character(coastal.flow.sites$Station[i])
  coastal.q=rbind(coastal.q,tmp)
  print(i)
}
coastal.q=merge(coastal.q,coastal.flow.sites,"Station")
coastal.q$Date.EST=date.fun(coastal.q$Date)
coastal.q2=merge(coastal.q,analysis.periods,by.x="Date.EST",by.y="DATE",all.x=T)

library(hydrostats)
hydrostat(subset(coastal.q2,Rivername=="Shark")$Data.Value)



subset(coastal.q2,Storm=="Wilma"&Station=="HARNEYFLAM")

plot(Data.Value~Date,subset(coastal.q,Station=="SHARKRIVBG"),type="l",xlim=wilma.period,ylab="Discharge (cfs)",lwd=1.5)
with(subset(coastal.q,Station=="HARNEYFLAM"),lines(Date,Data.Value,col="red",lwd=1.5))
legend("topright",legend=c("Shark","Harney"),lty=1,col=c("black","red"),ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)


with(subset(coastal.q,Station=="02290878"),lines(Date,Data.Value,col="blue"))
abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")

coastal.q.xtab=cast(coastal.q2,Date.EST+Storm+Period~RiverName,value="Data.Value",mean)
coastal.q.xtab$WY=WY(coastal.q.xtab$Date.EST)
range(coastal.q.xtab$WY)
library(mblm)
library(gvlma)
library(car)

plot(Harney~Shark,coastal.q.xtab);abline(0,1)

harney.fill.lm=lm(Harney~Shark,coastal.q.xtab)
gvlma(harney.fill.lm)
layout(matrix(1:4,2,2));plot(harney.fill.lm)
dev.off()
shapiro.test(residuals(harney.fill.lm));hist(residuals(harney.fill.lm))

harney.fill.lm=lm(Harney~Shark,subset(coastal.q.xtab,is.na(Harney)==F&is.na(Shark)==F))
mod.boot=Boot(harney.fill.lm,R=100, method="residual")
mod.boot
hist(mod.boot)

wy.samp=sample(2005:2018,5)
plot(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp));abline(0,1)
test=lm(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp))
gvlma(test)

#bootstrap regression
set.seed(3244)
bstar=NULL
n=nrow(coastal.q.xtab)
B=200
for(draw in 1:B){
  Dstar=subset(coastal.q.xtab,WY%in%sample(2005:2018,5))
  mod=lm(Harney~Shark,Dstar)
  bstar=rbind(bstar,coef(mod))
}
bstar

bstar=NULL
wys=2005:2018
n=length(wys)
for(i in 1:n){
  Dstar=subset(coastal.q.xtab,WY==wys[i])
  mod=lm(Harney~Shark,Dstar)
  bstar=rbind(bstar,coef(mod))
}
bstar


plot(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp),type="n")
for(i in 1:nrow(bstar)){
  abline(a=bstar[i,1],b=bstar[i,2],col=adjustcolor("grey",0.5))
}

##
coastal.q.xtab$Harney.fill=with(coastal.q.xtab,ifelse(is.na(Harney),(harney.fill.lm$coefficients[2]*Shark)+harney.fill.lm$coefficients[1],Harney))


coastal.q.xtab[is.na(coastal.q.xtab)]=0
#coastal.q.xtab$Harney=with(coastal.q.xtab,ifelse(Harney==0,Shark,Harney))

plot(Shark~Date.EST,coastal.q.xtab,type="l",xlim=wilma.period,ylim=c(-1500,10000))
with(coastal.q.xtab,shaded.range(Date.EST,rep(-2000,length(Date.EST)),Shark,"grey"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark,Shark+Harney.fill,"blue"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark+Harney.fill,Shark+Harney.fill+Broad,"red"))

plot(Shark~Date.EST,coastal.q.xtab,type="l",xlim=irma.period,ylim=c(-1500,20000))
with(coastal.q.xtab,shaded.range(Date.EST,rep(-2000,length(Date.EST)),Shark,"grey"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark,Shark+Harney,"blue"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark+Harney,Shark+Harney+Broad,"red"))



plot(Data.Value~Date,subset(coastal.q,Station=="SHARKRIVBG"),type="l",xlim=irma.period)
with(subset(coastal.q,Station=="HARNEYFLAM"),lines(Date,Data.Value,col="red"))
with(subset(coastal.q,Station=="02290878"),lines(Date,Data.Value,col="blue"))
abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")

plot(Data.Value~Date,subset(coastal.q,Station=="TAYLORS3"),type="l",xlim=irma.period)
with(subset(coastal.q,Station=="EASTCK"),lines(Date,Data.Value,col="red"))
with(subset(coastal.q,Station=="MUD_CRKM"),lines(Date,Data.Value,col="blue"))

abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")





### 
##
##
##
##
##
##
##

coastal.q2=coastal.q
coastal.q2$Q=coastal.q2$Data.Value
coastal.q2$Q=with(coastal.q2,ifelse(Data.Value<0,0,Data.Value))

tmp.dat=ts.format(subset(coastal.q2,Storm=="Wilma"&RiverName=="Shark")[,c("Date","Q")],format="%Y-%m-%d")
tmp.dat=ts.format(subset(coastal.q2,RiverName=="Shark")[,c("Date","Q")],format="%Y-%m-%d")
high.spells.v2(tmp.dat,hydro.year=T)
quantile(tmp.dat$Q,probs=c(0.9),na.rm=T)



##
##
##
##



subset(coastal.q2,Storm=="Wilma"&Station=="HARNEYFLAM")

plot(Data.Value~Date,subset(coastal.q,Station=="SHARKRIVBG"),type="l",xlim=wilma.period,ylab="Discharge (cfs)",lwd=1.5)
with(subset(coastal.q,Station=="HARNEYFLAM"),lines(Date,Data.Value,col="red",lwd=1.5))
legend("topright",legend=c("Shark","Harney"),lty=1,col=c("black","red"),ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)


with(subset(coastal.q,Station=="02290878"),lines(Date,Data.Value,col="blue"))
abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")

coastal.q.xtab=cast(coastal.q2,Date.EST+Storm+Period~RiverName,value="Data.Value",mean)
coastal.q.xtab$WY=WY(coastal.q.xtab$Date.EST)
range(coastal.q.xtab$WY)
library(mblm)
library(gvlma)
library(car)

plot(Harney~Shark,coastal.q.xtab);abline(0,1)

harney.fill.lm=lm(Harney~Shark,coastal.q.xtab)
gvlma(harney.fill.lm)
layout(matrix(1:4,2,2));plot(harney.fill.lm)
dev.off()
shapiro.test(residuals(harney.fill.lm));hist(residuals(harney.fill.lm))

harney.fill.lm=lm(Harney~Shark,subset(coastal.q.xtab,is.na(Harney)==F&is.na(Shark)==F))
mod.boot=Boot(harney.fill.lm,R=100, method="residual")
mod.boot
hist(mod.boot)

wy.samp=sample(2005:2018,5)
plot(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp));abline(0,1)
test=lm(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp))
gvlma(test)

#bootstrap regression
set.seed(3244)
bstar=NULL
n=nrow(coastal.q.xtab)
B=200
for(draw in 1:B){
  Dstar=subset(coastal.q.xtab,WY%in%sample(2005:2018,5))
  mod=lm(Harney~Shark,Dstar)
  bstar=rbind(bstar,coef(mod))
}
bstar

bstar=NULL
wys=2005:2018
n=length(wys)
for(i in 1:n){
  Dstar=subset(coastal.q.xtab,WY==wys[i])
  mod=lm(Harney~Shark,Dstar)
  bstar=rbind(bstar,coef(mod))
}
bstar


plot(Harney~Shark,subset(coastal.q.xtab,WY%in%wy.samp),type="n")
for(i in 1:nrow(bstar)){
  abline(a=bstar[i,1],b=bstar[i,2],col=adjustcolor("grey",0.5))
}

##
coastal.q.xtab$Harney.fill=with(coastal.q.xtab,ifelse(is.na(Harney),(harney.fill.lm$coefficients[2]*Shark)+harney.fill.lm$coefficients[1],Harney))


coastal.q.xtab[is.na(coastal.q.xtab)]=0
#coastal.q.xtab$Harney=with(coastal.q.xtab,ifelse(Harney==0,Shark,Harney))

plot(Shark~Date.EST,coastal.q.xtab,type="l",xlim=wilma.period,ylim=c(-1500,10000))
with(coastal.q.xtab,shaded.range(Date.EST,rep(-2000,length(Date.EST)),Shark,"grey"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark,Shark+Harney.fill,"blue"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark+Harney.fill,Shark+Harney.fill+Broad,"red"))

plot(Shark~Date.EST,coastal.q.xtab,type="l",xlim=irma.period,ylim=c(-1500,20000))
with(coastal.q.xtab,shaded.range(Date.EST,rep(-2000,length(Date.EST)),Shark,"grey"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark,Shark+Harney,"blue"))
with(coastal.q.xtab,shaded.range(Date.EST,Shark+Harney,Shark+Harney+Broad,"red"))



plot(Data.Value~Date,subset(coastal.q,Station=="SHARKRIVBG"),type="l",xlim=irma.period)
with(subset(coastal.q,Station=="HARNEYFLAM"),lines(Date,Data.Value,col="red"))
with(subset(coastal.q,Station=="02290878"),lines(Date,Data.Value,col="blue"))
abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")

plot(Data.Value~Date,subset(coastal.q,Station=="TAYLORS3"),type="l",xlim=irma.period)
with(subset(coastal.q,Station=="EASTCK"),lines(Date,Data.Value,col="red"))
with(subset(coastal.q,Station=="MUD_CRKM"),lines(Date,Data.Value,col="blue"))

abline(v=c(wilma.landfall.FL,irma.landfall.FL),col="red")

