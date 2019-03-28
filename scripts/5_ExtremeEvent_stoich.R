# Stoich ------------------------------------------------------------------

fce.grahl.auto=read.csv(paste0(data.path,"/WQ/LT_ND_Grahl_001.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.losada.auto=read.csv(paste0(data.path,"/WQ/LT_ND_Losada_001.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rondeau.auto=read.csv(paste0(data.path,"/WQ/LT_ND_Rondeau_001.txt"),na.strings=c("-9999","-9999.00","-9999.000"))
fce.rubio.auto=read.csv(paste0(data.path,"/WQ/LT_ND_Rubio_001.txt"),na.strings=c("-9999","-9999.00","-9999.000"))

names(fce.grahl.auto)
names(fce.losada.auto)
names(fce.rondeau.auto)
names(fce.rubio.auto)

fce.wq.auto=rbind(fce.grahl.auto,fce.losada.auto,fce.rondeau.auto,fce.rubio.auto)
rm(fce.grahl.auto,fce.losada.auto,fce.rondeau.auto,fce.rubio.auto)
fce.wq.auto$Date.EST=date.fun(as.character(fce.wq.auto$DATE),form="%Y-%m-%d")
fce.wq.auto$DATE=date.fun(as.character(fce.wq.auto$DATE),form="%Y-%m-%d")

unique(fce.wq.auto$SITENAME)
fce.wq.sites=data.frame(SITENAME=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Station.ID=c(paste0("SRS",c("1a","1c","1d",2,3,4,5,6)),paste0("TS/PH",c("1a","1b",2,3,"6b","6a","7a","7b"))),
                        Region=c(rep("SRS",8),rep("TS",8)))
fce.wq.auto=merge(fce.wq.auto,fce.wq.sites,"SITENAME")
unique(fce.wq.auto$SITENAME)

fce.wq.auto$TP.ugL=round((fce.wq.auto$TP*P.mw),4);#convert uM concentration to ug/L
fce.wq.auto$TN.mgL=round((fce.wq.auto$TN*N.mw)*0.001,2);#convert uM concentration to mg/L
fce.wq.auto$log.TP=log(fce.wq.auto$TP)
fce.wq.auto$log.TN=log(fce.wq.auto$TN)

fce.wq.auto2=merge(fce.wq.auto,analysis.periods,"DATE",all.x=T)

srs.wilma.sma=sma(log.TN~log.TP*Period,subset(fce.wq.auto2,Region=="SRS"&Storm=="Wilma"),slope.test = 1)
plot(srs.wilma.sma)
plot(srs.wilma.sma,which="residual")
plot(srs.wilma.sma,which="qq")
summary(srs.wilma.sma)

srs.irma.sma=sma(log.TN~log.TP*Period,subset(fce.wq.auto2,Region=="SRS"&Storm=="Irma"),slope.test = 1)
plot(srs.irma.sma)
plot(srs.irma.sma,which="residual")
plot(srs.irma.sma,which="qq")
summary(srs.irma.sma)






## check FCE_WQ_V1.r for begining SMA analysis

sites=c(paste0("SRS",c("1d",2,3,4,5,6)),paste0("TS/PH",c("1a",2,3,"6a","7a")))
SRS.site.sma=data.frame()
SRS.site.sma.compare=data.frame()
storm.val=c("Wilma","Irma")
# SRS Only
#storm (Wilma only)
#j=1;
for( j in 1:2){
  #tiff(filename=paste0(plot.path,"/SRS_",storm.val[j],"_sma.tiff"),width=7,height=7,units="in",res=220,type="windows",compression=c("lzw"),bg="white")
  par(family="serif",mar=c(1.5,2,1.5,0.5),oma=c(3,1.5,0.5,0.25))
  layout(matrix(1:18,6,3,byrow=T))
  #station
  for(i in 2:6){
    site.sma=sma(log.TN~log.TP*Period,subset(fce.wq.auto2,Station.ID==sites[i]&Storm==storm.val[j]),slope.test = 1)
    
    rslts=data.frame(Storm=storm.val[j],SITE=sites[i],Periods=site.sma$groups,
                     r2.val=as.numeric(unlist(site.sma$r2)),
                     slope=with(site.sma$coef,c(Post[2,1],Pre[2,1])),
                     inter=with(site.sma$coef,c(Post[1,1],Pre[1,1])),
                     F.val=unlist(site.sma$slopetest)[c(1,8)],
                     r.val=unlist(site.sma$slopetest)[c(2,9)],
                     p.val=unlist(site.sma$slopetest)[c(3,10)])
    rslts.slopecomp=data.frame(Storm=storm.val[j],SITE=sites[i],
                               LR=site.sma$commoncoef$LR,
                               p.val=site.sma$commoncoef$p)
    SRS.site.sma=rbind(SRS.site.sma,rslts)
    SRS.site.sma.compare=rbind(SRS.site.sma.compare,rslts.slopecomp)
    
    plot(site.sma)
    plot(site.sma,which="residual");mtext(side=3,sites[i])
    plot(site.sma,which="qq")
  }
}
