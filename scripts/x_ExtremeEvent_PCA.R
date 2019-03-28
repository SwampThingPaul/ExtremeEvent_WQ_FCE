
# PCA ---------------------------------------------------------------------

vars=c("Date","WQSite","Data.Value.m","depth.cm","hydroperiod.90dfreq")
wq.wl=merge(grab.wq2,WL.dat[,vars],by.x=c("Station.ID","Date.EST"),by.y=c("WQSite","Date"))
wq.wl=merge(wq.wl,WL.dat.xtab.melt,by.x=c("Date.EST","Region"),by.y=c("Date","Region"))
wq.wl=merge(wq.wl,rf.dat2[,c("Date.EST","TRF.m","cum.rf")],"Date.EST")

plot(TP~cum.rf,subset(wq.wl,Region=="SRS"))
names(wq.wl)

## PCA Analysis
my.rda=rda(tmp.dat[,-1:-4],scale=T)
biplot(my.rda)

eig <- my.rda$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)

#tiff(filename=paste0(plot.path,"/PCA_screplots.tiff"),width=6.75,height=3,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(0.5,1.5,1,1),oma=c(3.5,2,0.25,0.25),mgp=c(3,1,0));
layout(matrix(seq(1,2,1),1,2,byrow=T));

ylim.val=c(0,105);by.y=20;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)
points(x,eig.pca.all$cumvariance,pch=21,bg="indianred1",cex=1.25)
axis_fun(1,line=-0.7,x,x,seq(1,9,1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=2,"Percentage of Variances")
leg.x=x[1,]+(x[length(vars)-1,]-x[1,])/2
leg.y=ylim.val[1]-25
legend.text=c("Absolute","Cumulative")
pt.col=c("grey","indianred1")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1.55,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

ylim.val=c(0,3);by.y=1;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2)
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n")
abline(h=ymaj,lty=3,col="grey")
x=barplot(eig.pca$eig,ylim=ylim.val,col="grey",yaxt="n",add=T)
axis_fun(1,line=-0.7,x,x,seq(1,9,1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Principal Components")
mtext(side=2,line=1.5,"Eigenvalue")
dev.off()

plot(c(-2,2),c(-2,2),type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
x=ordiellipse(my.rda,group=tmp.dat$Region,draw="polygon",label=F,col=c("dodgerblue1","indianred1"),border=T,cex=0.8)

labs.rownames=merge(data.frame(vars=rownames(var.coord)),rownames.xwalk,sort=F)$short
scrs=scores(my.rda,display=c("sites","species"))

plot(c(-2,2),c(-2,2),type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
x=ordiellipse(my.rda,group=tmp.dat$Region,draw="polygon",label=F,col=c("indianred1","dodgerblue1"),border=T,cex=0.8)

labs.val=c("SRS","TS")
ellipse.center=data.frame()
for(i in 1:length(labs.val)){
  tmp=data.frame(Region=labs.val[i],pc1=as.numeric(x[[labs.val[i]]]$center[1]),pc2=as.numeric(x[[labs.val[i]]]$center[2]))
  ellipse.center=rbind(ellipse.center,tmp)
}
label_line=function(cen.x,cen.y,x2,y2,label,lty=1,cex=1,col="black",lwd=1,font=1,pos=NULL,pt.plot=T){
  if(pt.plot==F){NA}else{points(cen.x,cen.y,pch=19,col=col,cex=cex)}
  lines(c(cen.x,x2-(x2*0.06)),c(cen.y,y2-(y2*0.06)),lty=lty,col=col,lwd=lwd)
  text(x2,y2,label=label,cex=cex,col=col,font=font,pos=pos)
}

#tiff(filename=paste0(plot.path,"PCA_FCE.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.5,1.75),oma=c(2,1.75,0.25,0),mgp=c(3,1,0));
#layout(matrix(seq(1,2,1),1,2,byrow=T));

xlim.val=c(-3,3);by.x=1;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-3,3);by.y=1;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA)
abline(h=0,v=0,lty=3,col="grey")
with(tmp.dat,points(scrs$sites,pch=ifelse(tmp.dat$Region=="SRS",21,22),bg=adjustcolor(ifelse(tmp.dat$Region=="SRS","grey","white"),0.5),col=adjustcolor(ifelse(tmp.dat$Region=="SRS","grey30","black"),0.5),cex=1,lwd=0.5))
arrows(0,0,scrs$species[,1],scrs$species[,2],length = 0.1, angle = 15, code = 2,col="indianred1",lwd=1.5)
with(scrs,text(species[,1]+0.25,species[,2],labels=labs.rownames,cex=0.75))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1)
axis_fun(2,ymaj,ymin,format(ymaj),1)
mtext(side=1,line=1.8,paste0("PCA 1 (",round(eig.pca$variance[1],1),"%)"))
mtext(side=2,line=2,paste0("PCA 2 (",round(eig.pca$variance[2],1),"%)"))
legend("topleft",legend=c("SRS","TS"),pch=c(21,22),lty=0,pt.bg=adjustcolor(c("grey","white"),0.5),col=adjustcolor(c("grey30","black"),0.5),pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()




