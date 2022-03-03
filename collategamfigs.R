

#fig arrange

#run TPfull, TPfullligh, TPfullTP


mscale=c(0,
         max(pdat$upper, na.rm=T)+sscal)
#####

###tiff('baseR_figure.tiff', width =130, 
##  height = 80, pointsize = 12,
## units = 'mm', res = 300)

###jpeg('baseR_figure.jpeg', width =130,height = 80,
## units = 'mm', pointsize = 12,
##quality = 75,
##res = 300)

panlab=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(j)", "(i)" )


par(mfrow=c(3,3), mar=c(0.5,6,1,1), oma=c(4,1,3,2),
    cex.axis=1.1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
mmlab=c("TP", "kd", "chla")
gmlist=list(ft,fl,fc)
datlist=list(tpdat, lpdat,cpdat)
j=3
for (i in 1:3) {

pdat=as.data.frame(datlist[[i]])
fmod=gmlist[[i]]
yearfit=predict(fmod, newdata=pdat, se.fit=TRUE)
ffit=fitted_values(fmod)
pdat=cbind(pdat,ffit)


 # ylim1=c(0.5,2.25)
#  if (i==3) ylim1=ylim1/2
  datre=pdat[pdat$Station_Acronym==stvec[j],]
  plot(metric~year, xaxt="n",pch=16, data=datre,bty='l',cex=.8,
       ann=FALSE,  xlim=c(1972.5,2007.5),
       #ylim=ylim1,
       las=2)
  lines(fitted~year, lwd=2,data=datre,col=1)
  # lines(upper~year, data=datre,col="purple")
  #  lines(lower~year, data=datre,col="purple")
  x=c(rev(datre$lower), datre$upper)
  y=c(rev(datre$year), datre$year)
  polygon(y,x, col=adjustcolor("grey",alpha.f=0.3), border=NA )
  legend("topright",panlab[i], bty="n",
         inset=c(0.07,0), 
         cex=1.2)
  #  text(y=35, x=1990,stlabel[i], bty="n",
  #      cex=1.25 )
  # title(stlabel[i], cex=0.8)
 # mtext(stlabel[i], side = 3, line = -1, outer = FALSE, at = NA,
#        adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
  if ( i==3)axis(side=1, las=1)
  
mtext(mmlab[i], side = 2, line = 3.5, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
if (i==3) mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)

pout=partial_residuals(fmod, select = "s(year)", partial_match=TRUE)
sout=smooth_estimates(fmod, smooth = "s(year)", partial_match=TRUE)
pout=pout[,c(1,3,2)]

colnames(pout)=c("s(year):B", "s(year):HB", "s(year):C")

pdat=cbind(pdat,pout)
wind=which(pdat$Station_Acronym==stvec[j])
datsm=sout[sout$Station_Acronym==stvec[j],]
datsres=pdat[wind,]
yname=(paste0("s(year):",stvec[j]))
plot(datsres[,yname]~year, xaxt="n",pch=16,col="grey", 
     data=datsres,bty='l', cex=.8,
     ann=FALSE, xlim=c(1972.5,2007.5),
     #ylim=ylim2,
     las=2)
lines(est~year, lwd=2,data=datsm,col="dark green")
datsm$upper = datsm$est + (dfsm[i,2] * datsm$se)
datsm$lower = datsm$est - (dfsm[i,2] * datsm$se)
datsm$upper =  datsm$est + 2*datsm$se
datsm$lower =  datsm$est - 2*datsm$se
lines(upper~year, data=datsm,col="dark green")
lines(lower~year, data=datsm,col="dark green")
x=c(rev(datsm$lower), datsm$upper)
y=c(rev(datsm$year), datsm$year)
polygon(y,x, col=adjustcolor("dark green",alpha.f=0.1), border=NA )

if (i==2) mtext("smooth", side = 2, line = +3, outer = FALSE, at = NA,
                adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)

if ( i==3)axis(side=1, las=1)
legend("topright",panlab[i+3], bty="n",
       inset=c(0.07,0), 
       cex=1.2)
if (i==1) mtext(stlabel[j], side = 3, line = 0, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)

dtall=as.data.frame(derivatives(fmod, interval="simultaneous"))
dtall$smooth=as.factor(dtall$smooth)
dlit=split(dtall, dtall$smooth)


dseq=c(2,4,3)
if (i==1) dat=as.data.frame(dlit[dseq[j]-1])
if (i!=1)   dat=as.data.frame(dlit[dseq[j]])
#yaxl=c(min(dat$derivative), max(dat$derivative))
#if (yaxl[2]<0) yaxl[2]=0.01
yaxl=c(-0.001,0.0001)

if (i>1) yaxl=yaxl*5
colnames(dat)=colnames(dtall)

plot(derivative~data, data=dat,xaxt="n",ann=FALSE,bty='l',las=2,
 # ylim=c(-0.001,0.0001),
     typ="l", lwd=2, col="dark red", xlim=c(1972.5,2007.5))
lines(lower~data, data=dat, col="dark red")
lines(upper~data, data=dat,col="dark red")
x=c(rev(dat$lower), dat$upper)
#x[x>yaxl[2]]<-yaxl[2]+.01
#x[x< yaxl[1]]<-yaxl[1]-.01
y=c(rev(dat$data), dat$data)
polygon(y,x, col=adjustcolor("dark red",alpha.f=0.1), border=NA )

abline(h=0, col="dark grey", lwd=2)
abline(v=1994, col="dark grey", lwd=2)
if (i==2) mtext("derivative of smooth", side = 2, line = +4, outer = FALSE, at = NA,
                adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)

if ( i==3)axis(side=1, las=1)
legend("topright",panlab[i+6], bty="n",
       inset=c(-0.04,0), 
       cex=1.2)

}