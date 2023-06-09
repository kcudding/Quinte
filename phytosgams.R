rm(list=ls())

library("readxl")
library(car)
#install.packages("mgcv")
#this package fits a GAM
library("mgcv")
#install.packages("devtools") 
#library("devtools")
#devtools::install_github("gavinsimpson/gratia")
#this package calcuates the finite difference of a GAM curve
library("gratia")

#read in data supplied by DFO, rename and aggregate into annuual means
library(openxlsx)

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



form1="log(metric) ~  year+s(TP,bs='tp', m=2)"
#form4="log(metric) ~  t2(TP,year, bs=c('tp','tp'),m=2, full=TRUE)" #s
#form5="metric ~  te(TP,year, bs=c('tp', 'tp'), m=2)" #i

form2="(metric) ~  s(TP,bs='tp',m=2)+s(year,bs='tp',m=2)"
form3="(metric) ~  TP+s(year,bs='tp',m=2)"
form4="(metric) ~  s(TP)"
form5="(metric) ~  s(year)"
fcvec=c(form1,
        form2,  form3, form4, form5)

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
  
  
  agg <- aggregate(list(metric1=my_data[mmetric],TP=my_data$TP), 
                   by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                   FUN=mean, na.rm=TRUE, na.action="na.pass")
  agg2 <- aggregate(list(metricN=my_data[mmetric]), 
                    by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                    FUN=length)
  agg$n=agg2[,3]
  agg$weights <- agg$n/mean(agg$n)
  colnames(agg)[3]="metric1"
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
 #   xyear[1]=1984
  #  xyear[2]=2015
  cdat=cdatc
  cdat$metric=cdat[,metric2]
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
  
  N=1000 #number of points at which to evaluate the smooth
  
  rm(mm)
  mm=list()
  sc=vector()
  swf=vector()
  pdat=merge(cdat,spydat)
  
  lfm=lm(log(metric)~TP, data=pdat, weights=weights)
  lfm2=gam(log(metric)~s(TP) ,
      data = pdat, weights=weights, 
      
      na.action=na.omit,
      niterPQL=100,
      #correlation=corAR1(form = ~ year|Station_Acronym),
      method = "REML")
 # print(metric2)
#  print(summary(lfm))
  
  
  for (fr in 1:length(fcvec)){
   # print(fr)
    m<- gam(as.formula(fcvec[fr]) ,
            data = pdat, weights=weights, 
            
            na.action=na.omit,
            niterPQL=100,
            #correlation=corAR1(form = ~ year|Station_Acronym),
            method = "REML")
    sc[fr]=AIC(m)
    sw=shapiro.test(resid(eval(m)))
    swf[fr]=sw$p.value
    
  }
  print(sc-min(sc))
  print(swf)
#row.names(sc)=fcvec
  ind=which.min(sc)
  print(paste(ind, metric2))
 ind=5
 pdat$metric=(lfm2$residuals)
  fc<- gam(as.formula(fcvec[ind]) ,
           data = pdat, weights=weights, 
           
           na.action=na.omit,
           niterPQL=100,
           #correlation=corAR1(form = ~ year|Station_Acronym),
           
           method = "REML")
  cpdat=pdat
  summary(fc)
  AIC(fc)
  concurvity(fc, full=TRUE)
  print(round(concurvity(fc, full=FALSE)$worst,2))
  
  
  
  draw(fc, residuals=TRUE, rug=FALSE)
 # drt=plot(fc, seWithMean=TRUE, shift=coef(fc)[1], rug=FALSE, 
  #         +          shade=TRUE,residuals=TRUE, pch=16)
  fcfit=fitted_values(fc)
  pdat=cbind(pdat,fcfit)
  print(durbinWatsonTest(fc$residuals))
  fcderiv=as.data.frame(derivatives(fc, interval="simultaneous"))
  fcderiv$smooth=as.factor(fcderiv$smooth)
  lfc=split(fcderiv, fcderiv$smooth)
  par(mfrow=c(3,1), mar=c(1,4,1,3), oma=c(5,4,2,2),
      cex.axis=0.8, cex.lab=1.75, cex.main=1.2, cex.sub=1)
  dseq=c(1)
  for (i in seq_along(dseq)) {
    dat=as.data.frame(lfc["s(year)"])
    colnames(dat)=colnames(fcderiv)
    ymax=max(dat[, "upper"])
    ymin=min(dat[, "lower"])
    plot(fc, seWithMean=TRUE, shift=coef(fc)[1], rug=FALSE, 
         shade=TRUE,residuals=TRUE, pch=16, main=metric2)
    legend("right", legend=c(round(summary(fc)$r.sq,2),
                             round(summary(fc)$s.pv,2), 
                             round(summary(lfm)$adj.r.squared,2),
                             round(summary(lfm2)$r.sq,2)),
           bty="n", 
           lty=c(1,2,3),col=c(1,2,3) )
    abline(v=1994)
    plot(derivative~data, data=dat,
         ylim=c(ymin,ymax), 
         typ="l", lwd=2)
    lines(lower~data, data=dat)
    lines(upper~data, data=dat)
    abline(h=0)
    abline(v=1994)
  #  acf(fc$residuals)
    plot(lfm2, seWithMean=TRUE, shift=coef(fc)[1], rug=FALSE, 
         shade=TRUE,residuals=TRUE, pch=16, main=metric2)
    
    
  }
  
  }
