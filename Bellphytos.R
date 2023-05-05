library(strucchange)
library(ecotraj)

rm(list=ls())
#read in data supplied by DFO, rename and aggregate into annuual means
library(openxlsx)
library("readxl")
my_data <- read.xlsx("F:/homeoffice/kim/Quinte/QUINTE.xlsx", sheet = 1)

#rename col in recent file
colnames(my_data)[colnames(my_data)=="Epar#Ship"] <- "epar"
colnames(my_data)[colnames(my_data)=="Secchi_depth#Ship"] <- "secchi"

phytos=c("Unknown.algae.biomass","Cyanophyta.biomass",  "Chlorophyta.biomass",
         "Euglenophyto.biomass", "Chrysophyceae.biomass","Diatomeae.biomass" ,
         "Cryptophyta.biomass" , "Dinophyceae.biomass" )
#[65] "Cyanophyta.biomass"                             
#[66] "Chlorophyta.biomass"                            
#[67] "Euglenophyto.biomass"                           
#[68] "Chrysophyceae.biomass"                          
#[69] "Diatomeae.biomass"                              
#[70] "Cryptophyta.biomass"                            
#[71] "Dinophyceae.biomass" 


#mmetric="Chlorophyta.biomass" #select which metric to analyze

i <- c(64:71)  
my_data[ ,i] <- apply(my_data[ ,i], 2,function(x) as.numeric(as.character(x)))
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (ik in seq_along(phytos)) {
  mmetric=phytos[ik]
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
  
  xyear=c(1972,2008)
  spydat=spydat[spydat$year%in%c(1972:2015),]
  spydat$Station_Acronym=as.factor(spydat$Station_Acronym)

  site=stvec
  bxlab=c("","", "Breaks")
  bylab=c("","BIC", "")
  ly=expression(paste("Light Attentuated (", m^-1, ")"))
  lylab=c("", ly, "")
  
  
  #TP model
  tp1=(lm(metric~TP, weights=weights, data=spydat))
  tp2=(lm(metric~year, weights=weights, data=spydat))
  tp4=(lm(metric~TP+year, weights=weights, data=spydat))
  
 
  
  cop2=AIC(tp1,tp2,tp4)
  print(cop2)
  print(c(summary(tp1)$adj.r.squared,
          summary(tp2)$adj.r.squared,summary(tp4)$adj.r.squared))
  
  m1=which.min(cop2$AIC)
  mod=row.names(cop2)[m1]
  bestmod=as.name(mod)
  
  #summary(eval(bestmod))
  mlab="phyto"
  
  mpred=as.data.frame(predict(tp4, spydat, interval="confidence"))
  
  spydat=cbind(spydat, mpred)
  sscal=1
  #spydat$tres=spydat$TP-spydat$tmodel
  
  
    pdat=spydat
    
    mscale=c(min(pdat$metric, na.rm=TRUE)-sscal,
             max(pdat$metric, na.rm=T)+sscal)
    # mscale=c(0,4500)
    amscale=mscale
    mss=which(!complete.cases(pdat))
    if (length(mss)>0) pdat=pdat[-mss,]
  
      
      amscale[2]=mscale[2]/2
   
    ymax=round(max(pdat$metric, na.rm=T),2)
    plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
         ann=FALSE, xlim=xyear,
       #  ylim=amscale, 
      #   yaxp=c(0,ymax,3),
         las=2)
    
   
      axis(side=1, las=1)
      mtext(mmetric,side=1)
  
    
    legend("bottomleft",stlabel[i], bty="n",
           cex=1.5)
    
    lines(pdat$fit~pdat$year, col = 2, lwd=2)
    lines(pdat$lwr~pdat$year, col = "pink", lwd=2, lty=2)
    lines(pdat$upr~pdat$year, col = "pink", lwd=2, lty=2)
   
  
  
  treg=summary(tp4)$coefficients
  
} 
