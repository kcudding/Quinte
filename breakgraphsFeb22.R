##########BREAKPOINT

#read in data supplied by DFO, rename and aggregate into annuual means
library(openxlsx)

my_data <-(as.data.frame(read.csv("E:/homeoffice/kim/Quinte/Quinte2022.csv")))
#rename col in recent file
colnames(my_data)[1]<-"Station_Acronym"
colnames(my_data)[colnames(my_data)=="Kd"] <- "epar"
colnames(my_data)[colnames(my_data)=="Secchi"] <- "secchi"
colnames(my_data)[colnames(my_data)=="chl"] <- "Chla"

colnames(my_data)[colnames(my_data)=="TotalPhyto"]<-"Phyto_BM"

mmetric="Chla" #select which metric to analyze

#aggregate and run analysis on selected sites
frm=paste(mmetric,"~ Station_Acronym+year")
spydat=aggregate(as.formula(frm), my_data, length)
mdat=aggregate(as.formula(frm), my_data, mean)
colnames(spydat)[3]="n"
spydat[4]=mdat[3]
colnames(spydat)[4]="metric"
spydat$weights <- spydat$n/mean(spydat$n)


agg <- aggregate(list(metric=my_data[mmetric],TP=my_data$TP), 
                 by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                 FUN=mean, na.rm=TRUE, na.action="na.pass")
agg2 <- aggregate(list(metricN=my_data[mmetric]), 
                  by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                  FUN=length)
agg$n=agg2[,3]
agg$weights <- agg$n/mean(agg$n)
colnames(agg)[3]="metric"
colnames(agg)[4]="TP"

#select sites on which to run analysis
stvec=c("B", "HB",  "C")

stlabel=c("Belleville", "Hay Bay", "Conway")
#select sites on which to run analysis

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)


spydat=agg[agg$Station_Acronym%in%stvec,]

xyear=c(1972,2013)
spydat=spydat[spydat$year%in%c(1972:2013),]
spydat$Station_Acronym=as.factor(spydat$Station_Acronym)
str(spydat)

#rm(list=ls())
library(strucchange)




site=stvec
bxlab=c("","", "Breaks")
bylab=c("","BIC", "")
ly=expression(paste("Light Attentuated (", m^-1, ")"))
lylab=c("", ly, "")
site[i]

#TP model
tp1=(lm(TP~year:Station_Acronym, weights=weights, data=spydat))
tp2=(lm(TP~year+Station_Acronym, weights=weights, data=spydat))
tp3=(lm(TP~year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))


AIC(tp1,tp2,tp3)
mpred=as.data.frame(predict(tp3, spydat, interval="confidence"))
spydat$tmodel=mpred$fit
spydat$tlow=mpred$lwr
spydat$tup=mpred$upr
spydat$tres=spydat$TP-spydat$tmodel
par(mfrow=c(3,1), mar=c(0,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:length(stvec)) {
  pdat=spydat[spydat$Station_Acronym==site[i],]
  if (i==3) {
  mss=which(is.na(pdat$TP))
  pdat=pdat[-mss]
  amscale[2]=mscale[2]/3
  }
ymax=round(max(pdat$TP, na.rm=T),2)
plot(TP~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
     ann=FALSE, xlim=xyear, ylim=c(0,ymax+0.01), las=2, yaxp=c(0,ymax,3))

if ( i==3)axis(side=1, las=1)
legend("bottomleft",stlabel[i], bty="n",
       cex=1.5)

lines(pdat$tmodel~pdat$year, col = 2, lwd=2)
lines(pdat$tlow~pdat$year, col = "pink", lwd=2, lty=2)
lines(pdat$tup~pdat$year, col = "pink", lwd=2, lty=2)
}

for (i in 1:length(stvec)) {
  pdat=spydat[spydat$Station_Acronym==stvec[i],]
  dat=ts(spydat$metric[spydat$Station_Acronym==site[i]], start=spydat$year[1], 
         end=spydat$year[length(spydat$year)])
  tt=1:length(dat)
  #b=breakpoints(dat~1)
  b=breakpoints(metric ~ year+tres, data=pdat)
  #b=breakpoints(metric ~ year, data=pdat)
  print(b)
  mss=which(is.na(dat))
  if (sum(!is.na(b$breakpoints))>0){
  
  # fit segmented model
  fac.ri <- breakfactor(b, label = "seg")
   if (length(mss)>0){
  bfac.ri=fac.ri[1:(mss-1)]
  bfac.ri=c(bfac.ri, NA)
  bfac.ri=c(bfac.ri, fac.ri[(mss):length(fac.ri)])
   }else{
     bfac.ri= fac.ri
  }
  
  bfac.ri=as.factor(bfac.ri)
  fm.ri <- lm(metric ~ year + tres+ bfac.ri, data=pdat)
  summary(fm.ri)
  }
  amscale=mscale
  amscale=c(0, amscale[2])
  if (i==3) amscale[2]=mscale[2]/3
  plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
       ann=FALSE, xlim=xyear, las=2)
  
  if ( i==3)axis(side=1, las=1)
  legend("bottomleft",stlabel[i], bty="n",
         cex=1.5)

  bfm=predict( fm.ri, newdata=pdat)
  
  lines(bfm~pdat$year, col = 2, lwd=2)

  
  
  #convert confidence intervals from julian day to dates
  #and draw
  
  #draw breakpoints & confidence intervals 
 
  ci=confint(b)$confint
  hgt=10
  if (i==3) hgt=4
  yn=rep(hgt,nrow(ci))+rnorm(5,1,1)
  s=1:nrow(ci)
  xn=pdat$year[ci[,1]]
  xx=pdat$year[ci[,3]]
  yx=yn
  segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
  abline(v=pdat$year[ci[,2]], lty=2, col="blue")
  }

  
  
  
 
  #abline(v=summary(b)$breakdates[b,], col="red")
  
 # text(breakdates(b)+3,max(dat)-0.05*max(dat), breakdates(b))
  text(breakdates(b)+3,30, breakdates(b))
  mtext(lylab[i], side = 2, line = -2, outer = TRUE, at = NA,
        adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)
  
#  lines(fitted(b, breaks = sum(!is.na(b$breakpoints))), 
 #       col = 4)
if (sum(!is.na(b$breakpoints))>0) {
    lines(confint(b, breaks = sum(!is.na(b$breakpoints))), col="green")
  }
  
  

##########END BREAKPOINT