##########BREAKPOINT

rm(list=ls())
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
#mlab=expression(Light~attenuation~(m^-1)) 
mlab=expression(Chlorophyll~a~(mu~g~L^-1))
#mlab=expression(Total~phosophorus~(mg~L^-1))


sscal=1.25 #clorophyll and round 1

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


#metric model
tp1=(lm(metric~year:Station_Acronym, weights=weights, data=spydat))
tp2=(lm(metric~year+Station_Acronym, weights=weights, data=spydat))
tp3=(lm(metric~year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
#metric model
tp4=(lm(metric~TP+year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
tp5=(lm(metric~TP+year+Station_Acronym, weights=weights, data=spydat))
tp6=(lm(metric~TP:year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
tp7=(lm(metric~TP:Station_Acronym+year+Station_Acronym, weights=weights, data=spydat))
tp8=(lm(metric~TP+year:Station_Acronym, weights=weights, data=spydat))
tp9=(lm(metric~TP+year+Station_Acronym, weights=weights, data=spydat))
tp10=(lm(metric~TP:year:Station_Acronym, weights=weights, data=spydat))
tp11=(lm(metric~TP:Station_Acronym+year, weights=weights, data=spydat))



AIC(tp1,tp2,tp3, tp4,tp5,tp6, tp7, tp8, tp9, tp10, tp11)
best=tp5
mpred=as.data.frame(predict(tp5, spydat, interval="confidence"))

spydat=cbind(spydat, mpred)

#spydat$tres=spydat$TP-spydat$tmodel
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:length(stvec)) {
  pdat=spydat[spydat$Station_Acronym==stvec[i],]
  
  mscale=c(min(pdat$metric, na.rm=TRUE)-sscal,
           max(pdat$metric, na.rm=T)+sscal)
  amscale=mscale
  if (i==3) {
  mss=which(is.na(pdat$TP))
  pdat=pdat[-mss]
  amscale[2]=mscale[2]/3
  }
ymax=round(max(pdat$metric, na.rm=T),2)
plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
     ann=FALSE, xlim=xyear, ylim=c(0,ymax+0.01), las=2, yaxp=c(0,ymax,3))

if ( i==3)axis(side=1, las=1)
legend("bottomleft",stlabel[i], bty="n",
       cex=1.5)

lines(pdat$fit~pdat$year, col = 2, lwd=2)
lines(pdat$lwr~pdat$year, col = "pink", lwd=2, lty=2)
lines(pdat$upr~pdat$year, col = "pink", lwd=2, lty=2)
}


tiff('baseR_figure.tiff', width =85, 
     height = 180, pointsize = 12,
     units = 'mm', res = 300)

jpeg('baseR_figure.jpeg', width =85, 
     height = 180,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

par(mfrow=c(3,1), mar=c(0.5,3,1,1), oma=c(4,3,0,0),
    cex.axis=1.2, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:length(stvec)) {
  pdat=spydat[spydat$Station_Acronym==stvec[i],]
  dat=ts(spydat$metric[spydat$Station_Acronym==site[i]], start=spydat$year[1], 
         end=spydat$year[length(spydat$year)])
  tt=1:length(dat)
  mss=which(!complete.cases(pdat))
  if (length(mss)>0) pdat=pdat[-mss,]
  
  #b=breakpoints(dat~1)
  b=breakpoints(metric ~ year+TP, data=pdat)
  #b=breakpoints(metric ~ year, data=pdat)
  print(b)

  if (sum(!is.na(b$breakpoints))>0){
  
  # fit segmented model
  fac.ri <- breakfactor(b, label = "seg")
#   if (length(mss)>0){
#  bfac.ri=fac.ri[1:(mss-1)]
#  bfac.ri=c(bfac.ri, NA)
#  bfac.ri=c(bfac.ri, fac.ri[(mss):length(fac.ri)])
#   }else{
#     bfac.ri= fac.ri
  
  
  time=as.factor(fac.ri)
  fm.ri <- lm(metric ~ time+TP+year, data=pdat, weights=weights)
  fm.riN <- lm(metric ~ TP+year, data=pdat,weights=weights)
  AIC(fm.ri, fm.riN)
  summary(fm.ri)
  } else {
    fm.ri <- lm(metric ~ TP+year, data=pdat,weights=weights)
  }
  amscale=mscale
  amscale=c(0, amscale[2])
  if (i==3) amscale[2]=mscale[2]/3
  plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
       ann=FALSE, xlim=xyear, ylim=c(0,45),las=2,yaxp=c(0,round(40,2),4))
  
  if ( i==3)axis(side=1, las=1)
  legend("topright",stlabel[i], bty="n",
         cex=1.5)

  bfm=predict( fm.ri, newdata=pdat)
  
  lines(bfm~pdat$year, col = 2, lwd=2)

  
  #convert confidence intervals from julian day to dates
  #and draw
  
  #draw breakpoints & confidence intervals 
 
  ci=confint(b)$confint
  ci=ci[summary(fm.ri)$coefficients[2:(1+length(b$breakpoints)),4]<0.05,]
  hgt=40
#  if (i==3) hgt=4
  tm=nrow(ci)
  if (length(tm)==0) tm=1
  yn=rep(hgt,tm)
  
  if (tm==1) {
    s=1
    xn=pdat$year[ci[1]]
    xx=pdat$year[ci[3]]
    yx=yn
    segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
    text(pdat$year[ci[2]]-2.5,44, pdat$year[ci[2]])
    abline(v=pdat$year[ci[2]], lty=2, col="blue")
  } else {
    yn=rep(hgt,tm)
    
  s=1:tm
  
  xn=pdat$year[ci[,1]]
  xx=pdat$year[ci[,3]]
  yx=yn
  segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
  abline(v=pdat$year[ci[,2]], lty=2, col="blue")
  text(pdat$year[ci[,2]]-2.5,44, pdat$year[ci[,2]])
  }

  
  
  mtext(mlab, side = 2, line = 0, outer = TRUE, at = NA,
        adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)
  mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
        adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)
  

}
legend("right", c("data", "model", "breakpoint",
                  "breakpoint CI"),
       lty=c(1,1,2,1), col=c("grey", "red","blue", "blue"), 
       pch=c(16,NA,NA,NA),
       bty="n", cex=1.2)
dev.off()

##########END BREAKPOINT