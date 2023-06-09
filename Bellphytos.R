library(strucchange)
library(ecotraj)
library(car)
library(nlme)

rm(list=ls())
#read in data supplied by DFO, rename and aggregate into annuual means
library(openxlsx)
library("readxl")
my_data <- read.xlsx("F:/homeoffice/kim/Quinte/QUINTE.xlsx", sheet = 1)
my_data2 <- read.xlsx("F:/homeoffice/kim/Quinte/All phyto annual average Comp and Ind date avgwith count.xlsx", 
                     sheet=3, startRow=4)
my_data3 <- read.xlsx("F:/homeoffice/kim/Quinte/All phyto annual average Comp and Ind date avgwith count.xlsx", 
                      sheet=3, startRow=4)

my_data4 <- read.xlsx("F:/homeoffice/kim/Quinte/All phyto annual average Comp and Ind date avgwith count.xlsx", 
                      sheet=1, startRow=4)
#rename col in recent file
colnames(my_data)[colnames(my_data)=="Epar#Ship"] <- "epar"
colnames(my_data)[colnames(my_data)=="Secchi_depth#Ship"] <- "secchi"

phytos=c("Cyanophyta.biomass",  "Chlorophyta.biomass",
         "Chrysophyceae.biomass","Diatomeae.biomass" ,
         "Cryptophyta.biomass" , "Dinophyceae.biomass" , "Total.Phyto.biomass.searching.Excel.files")
phytos2=c("C.Cyano",   "C.Chloro", "C.Chryso",  "C.Diatom",  "C.Crypto", 
"C.Dino",    "C.Total")
phytos2I=c("I.Cyano",   "I.Chloro",  "I.Chryso",  "I.Diatom",  "I.Crypto",  "I.Dino",   
"I.Total")

#[65] "Cyanophyta.biomass"                             
#[66] "Chlorophyta.biomass"                            
#[67] "Euglenophyto.biomass"                           
#[68] "Chrysophyceae.biomass"                          
#[69] "Diatomeae.biomass"                              
#[70] "Cryptophyta.biomass"                            
#[71] "Dinophyceae.biomass" 


#mmetric="Chlorophyta.biomass" #select which metric to analyze

i <- c(64:72)  
my_data[ ,i] <- apply(my_data[ ,i], 2,function(x) as.numeric(as.character(x)))
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (ik in seq_along(phytos)) {
  mmetric=phytos[ik]
  metric2=phytos2[ik]
  metricI=phytos2I[ik]
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
  stvec=c("B")
  
  stlabel=c("Belleville")
  #select sites on which to run analysis
  
  
    spydat=agg[agg$Station_Acronym%in%stvec,]
    cdat=my_data2[my_data2$station%in%stvec,]
    cdatc=my_data4[my_data4$station%in%stvec,]
    
    
  xyear=c(max(min(cdat$year), min(spydat$year)),
          min(max(cdat$year), max(spydat$year)))
#  xyear[1]=1984
#  xyear[2]=2015
  cdat=cdatc
  spydat=spydat[spydat$year%in%c(xyear[1]:xyear[2]),]
  cdat=cdat[cdat$year%in%c(xyear[1]:xyear[2]),]
  spydat$Station_Acronym=as.factor(spydat$Station_Acronym)
  spydat[is.na(spydat) | spydat=="Inf"] = NA
  cdat[is.na(cdat) | cdat=="Inf"] = NA
  site=stvec
  bxlab=c("","", "Breaks")
  bylab=c("","BIC", "")
  ly=expression(paste("Light Attentuated (", m^-1, ")"))
  lylab=c("", ly, "")
  colnames(cdat)[35]="samples"
  
  #TP model
  tp1=(lm(metric~TP, weights=weights, data=spydat))
  tp2=(lm(metric~year, weights=weights, data=spydat))
  tp4=(lm(metric~TP+year, weights=weights, data=spydat))
  
  #TP model
  tp1=lm(log(cdat[,metric2])~spydat$TP, weights=spydat$weights, data=cdat)
  tp2=lm(log(cdat[,metric2])~spydat$year, weights=spydat$weights, data=cdat)
  tp4=lm(log(cdat[,metric2])~spydat$TP+spydat$year, weights=spydat$weights, data=cdat)
  
  cop2=AIC(tp1,tp2,tp4)
 
  print(mmetric)
  print(cop2)
  print(c(summary(tp1)$adj.r.squared,
          summary(tp2)$adj.r.squared,summary(tp4)$adj.r.squared))
  
  p1=pf(summary(tp1)$fstatistic[1],summary(tp1)$fstatistic[2],
     summary(tp1)$fstatistic[3],lower.tail=FALSE)
  p2=pf(summary(tp2)$fstatistic[1],summary(tp2)$fstatistic[2],
        summary(tp2)$fstatistic[3],lower.tail=FALSE)
  p4=pf(summary(tp4)$fstatistic[1],summary(tp4)$fstatistic[2],
        summary(tp4)$fstatistic[3],lower.tail=FALSE)
  print(c(p1,p2,p4))
  
  b=breakpoints(cdat[,metric2]~cdat$year)
  print(b)
  
  fac.ri <- breakfactor(b, label = "seg")
  #   if (length(mss)>0){
  #  bfac.ri=fac.ri[1:(mss-1)]
  #  bfac.ri=c(bfac.ri, NA)
  #  bfac.ri=c(bfac.ri, fac.ri[(mss):length(fac.ri)])
  #   }else{
  #     bfac.ri= fac.ri
  m1=which.min(cop2$AIC)
  mod=row.names(cop2)[m1]
  bestmod=as.name(mod)
  sw=shapiro.test(resid(eval(bestmod)))
  dw=durbinWatsonTest(eval(bestmod))
  
  if (sum(!is.na(b$breakpoints))>0) {
  time=as.factor(fac.ri)
 
  fm.ri <- lm(log(cdat[,metric2]) ~ time*year, data=cdat, weights=spydat$weights)
  fm.ri2 <- lm(log(cdat[,metric2]) ~ time+year, data=cdat, weights=spydat$weights)
 
    fm.ri3 <- lm(log(cdat[,metric2]) ~ time*spydat$TP, data=cdat, weights=spydat$weights)
    fm.ri4 <- lm(log(cdat[,metric2]) ~ time+spydat$TP, data=cdat, weights=spydat$weights)
  
    fm.ri5 <- lm(log(cdat[,metric2]) ~ time+year+spydat$TP, data=cdat, weights=spydat$weights)
    fm.ri6 <- lm(log(cdat[,metric2]) ~ time*year*spydat$TP, data=cdat, weights=spydat$weights)
    fm.ri7 <- lm(log(cdat[,metric2]) ~ time*year+spydat$TP, data=cdat, weights=spydat$weights)
    fm.ri8 <- lm(log(cdat[,metric2]) ~ time+year*spydat$TP, data=cdat, weights=spydat$weights)
    
  }
  bmods=print(AIC(fm.ri,fm.ri2,fm.ri3,fm.ri4,fm.ri5, fm.ri6, fm.ri7, fm.ri8))
  b1=which.min(bmods$AIC)
  bm1=row.names(bmods)[b1]
  bestmod=as.name(bm1)
  sw=shapiro.test(resid(eval(bestmod)))
  dw=durbinWatsonTest(eval(bestmod))
  print(summary(eval(bestmod)))
  


  
#  print(sw)
  dw=durbinWatsonTest(fm.ri)
  dw
  vif(eval(bestmod))
  mmet=cdat[,metric2]
  cmm=gls(log(mmet)~time*year,data=cdat, correlation = corAR1())
  
  cmm2=gls(log(mmet)~time+year,data=cdat, correlation = corAR1())
  v <- cmm$residuals
  attr(v,"std") <- NULL      # get rid of the additional attribute
  car::durbinWatsonTest( v )
  acf(v, main=mmetric)
  #summary(eval(bestmod))
  mlab="phyto"
  print(mmetric)
  print(AIC(fm.ri,cmm, cmm2))
  print(summary(cmm))
  mpred=as.data.frame(predict(eval(bestmod), 
              data.frame(year=cdat$year, TP=spydat$TP), interval="confidence"))
  mpred=as.data.frame(predict(fm.ri, 
      data.frame(year=cdat$year, time=time), interval="confidence"))
 # colnames(mpred)="fit"
  #mpred=exp(mpred)
  spydat=cbind(spydat, mpred)
  sscal=1
  #spydat$tres=spydat$TP-spydat$tmodel
  
  
    pdat=spydat
    
   # mscale=c(min(pdat$metric, na.rm=TRUE)-sscal,
    #         max(pdat$metric, na.rm=T)+sscal)
     mscale=c(0,4500)
 amscale=mscale
  #  mss=which(!complete.cases(pdat))
   # if (length(mss)>0) pdat=pdat[-mss,]
  
      
      amscale[2]=mscale[2]/2
   
    ymax=round(log(max(c(pdat$metric+0.001,cdat[,metric2]+0.001, cdat[,metricI]+0.001), na.rm=T)),2)
    ymin=round(log(min(c(pdat$metric+0.001,cdat[,metric2]+0.001, cdat[,metricI]+0.001),na.rm=T)),2)
    amscale[2]=ymax
    
    plot(log(cdat[,metric2])~cdat$year,xaxt="n",pch=16,col="grey",
         ann=FALSE, xlim=xyear,ylim=c(ymin,ymax+.1*ymax),
       #  ylim=amscale, 
      #   yaxp=c(0,ymax,3),
      type="b",
         las=2)
    title(main=mmetric)
    
   # text(cdat$year, cdat[,ik+2], text=cdat$samples, cex=.5)
    
    k=ik
    if (ik==7) k=8
   points(log(metric)~year, col="green", pch=15, type="b", data=spydat, cex=1.5)
   points(log(cdat[,metricI])~cdat$year, col="blue", pch=15, type="b", data=spydat)
     axis(side=1, las=1)
      mtext(mmetric,side=1)
  
    
    legend("bottomleft",stlabel[i], bty="n",
           cex=1.5)
    
    lines(mpred$fit~pdat$year, col = 2, lwd=2)
   # lines(pdat$lwr~pdat$year, col = "pink", lwd=2, lty=2)
  #  lines(pdat$upr~pdat$year, col = "pink", lwd=2, lty=2)
    text(paste(bestmod, round(summary(fm.ri)$adj.r.squared,2), 
               round(sw$p.value,2), dw$p), 
         x=2005, y=ymax-.1*ymax, cex=1)
    
  
  
  treg=summary(tp4)$coefficients

  
  #draw breakpoints & confidence intervals 
  if (sum(!is.na(b$breakpoints))>0){
    cis=confint(b)$confint
    cll=cis[summary(fm.ri)$coefficients[2:(1+length(b$breakpoints)),4]<1.,]
    ci=cis[summary(fm.ri)$coefficients[2:(1+length(b$breakpoints)),4]<0.05,]
    if (length((cll))!= 0) {
      if (is.null(nrow(cll))) abline(v=cdat$year[cll[2]], lty=2, col="light grey")
      else if (nrow(cll)>1) abline(v=cdat$year[cll[,2]], lty=2, col="light grey")
    }
    hgt=amscale[2]-(8/i)
    #  if (i==3) hgt=4
    tm=nrow(ci)
    if (length(tm)==0) {
      tm=1 
    } else if (tm==0)  {
      tm=1
    }
    yn=rep(hgt,tm)
    
    if (tm==1) {
      s=1
      xn=cdat$year[ci[1]]
      xx=cdat$year[ci[3]]
      yx=yn
      segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
      text(cdat$year[ci[2]]-2,amscale[2]-.5, cdat$year[ci[2]])
      abline(v=cdat$year[ci[2]], lty=2, col="blue")
      
      
    } else {
      yn=rep(hgt,tm)
      
      s=1:tm
      
      xn=cdat$year[ci[,1]]
      xx=cdat$year[ci[,3]]
      yx=yn
      segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
      abline(v=cdat$year[ci[,2]], lty=2, col="blue")
      text(cdat$year[ci[,2]]-2,amscale[2]-.5, cdat$year[ci[,2]])
     
    }
  }
  
 # if (i==1)legend("bottom", c("data", "model", "break",
  #                            "break CI", "ns break"), 
   #               lty=c(1,1,2,1,2), col=c("grey", "red","blue", "blue", "grey"), ncol=2,xjust=1, yjust=1,y.intersp = .75,
    #              pch=c(16,NA,NA,NA, NA),
     #             inset=c(0, -.025),
      #            bty="n", cex=1.2)
 # abline(v=1994)
  
}
