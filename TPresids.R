#install.packages("readxl")
rm(list=ls())
library("readxl")

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

mlab=expression(Light~attenuation~(m^-1)) 
#mlab=expression(Chlorophyll~a~(mu~g~L^-1))
#mlab=expression(Total~phosophorus~(mg~L^-1))
#rlab=expression(Rate~of~change~(mg~L^-1~year^-1))
rlab=expression(Rate~of~change~(m^-1~year^-1))
#rlab=expression(Rate~of~change~(mu~g~L^-1~year^-1))
sscal=0.0001 #phosophorus and round 3
sscal=1.25 #clorophyll and round 1
#sscal=0.005 #phosophorus and round 3
sscal=1
#aggregate and run analysis on selected sites
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
xyear=c(1978,2015)
spydat=agg[agg$Station_Acronym%in%stvec,]

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)
pdat=spydat
pdat=pdat[pdat$year%in%c(1978:2015),]
xyear=c(1978,2015)
str(pdat)
N=1000 #number of points at which to evaluate the smooth

pdat$Station_Acronym=as.factor(pdat$Station_Acronym)


par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)

str(pdat)
N=1000 #number of points at which to evaluate the smooth


form1="metric ~  Station_Acronym+s(TP,by=Station_Acronym, 
bs='tp',k=20,m=2)" 
form2="metric ~  s(TP,by=Station_Acronym, bs='tp',k=20,m=2)" 
form3="metric ~  s(TP,bs='tp',k=20,m=2)" 
form4="metric ~  Station_Acronym+s(TP,bs='tp',k=20,m=2)" 
frmvec1=c(form1,form2, form3, form4)
sc=vector()
for (fr in 1:length(frmvec1)){
  m<- gam(as.formula(frmvec1[fr]) ,
          data = pdat, weights=weights, 
          na.action=na.omit,
          niterPQL=100,
          #correlation=corAR1(form = ~ year|Station_Acronym),
          method = "REML")
  sc[fr]=AIC(m)
  
}

ind=which.min(sc)
ind=3
m<- gam(as.formula(frmvec1[ind]) ,
        data = pdat, weights=weights, 
        na.action=na.omit,
        niterPQL=100,
        #correlation=corAR1(form = ~ year|Station_Acronym),
        method = "REML")
summary(m)
AIC(m)
concurvity(m, full=TRUE)
round(concurvity(m, full=FALSE)$worst,2)

draw(m, residuals=TRUE, rug=FALSE)
TPfit=predict(m, newdata=pdat, se.fit=TRUE)
pdat$TPfit=TPfit$fit
pdat$TPfitse=TPfit$se
pdat$TPres=pdat$metric-pdat$TPfit
crit.t <- qt(0.975, df = df.residual(m)) 
pdat$TPupr = TPfit$fit + (crit.t *TPfit$se)
pdat$TPlwr = TPfit$fit + (crit.t *TPfit$se)

mscale=c(min(newYear$lower, na.rm=TRUE)-sscal,
         max(pdat$metric, na.rm=T)+sscal)




form6="TPres ~  Station_Acronym+s(year,by=Station_Acronym, bs='tp', k=20,m=1)" 
form7="TPres ~  s(year,by=Station_Acronym,bs='tp',k=20,m=2)"
form8="TPres ~  s(year,bs='tp',k=20,m=2)"




frmvec=c(form6,form7, form8)
rm(mm)
mm=list()
sc=vector()
for (fr in 1:length(frmvec)){
 fm<- gam(as.formula(frmvec[fr]) ,
          data = pdat, weights=weights, 
          
          na.action=na.omit,
          niterPQL=100,
          #correlation=corAR1(form = ~ year|Station_Acronym),
          method = "REML")
  sc[fr]=AIC(m)
  
  
  
}
ind=which.min(sc)
ind=2
fm<- gam(as.formula(frmvec[ind]) ,
         data = pdat, weights=weights, 
         
         na.action=na.omit,
         niterPQL=100,
         correlation=corAR1(form = ~ year|Station_Acronym),
         
         method = "REML")

summary(fm)
AIC(fm)
concurvity(fm, full=TRUE)
round(concurvity(fm, full=FALSE)$worst,2)
draw(m, residuals=TRUE)
draw(fm, residuals=TRUE, rug=FALSE)

plot(fm, seWithMean=TRUE, shift=coef(fm)[1], rug=FALSE, 
     shade=TRUE,residuals=TRUE, pch=16)
#by.resids=TRUE)

fmderiv=as.data.frame(derivatives(fm, interval="simultaneous"))
fmderiv$smooth=as.factor(fmderiv$smooth)
lfm=split(fmderiv, fmderiv$smooth)
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:3) {
  dat=as.data.frame(lfm[i])
  colnames(dat)=colnames(fmderiv)
 plot(derivative~data, data=dat, ylim=c(-1.1,1.1), typ="l", lwd=2)
 lines(lower~data, data=dat)
 lines(upper~data, data=dat)
 abline(h=0)
 abline(v=1994)
}


yearfit=predict(fm, newdata=pdat, se.fit=TRUE)
pdat$yrfit=yearfit$fit
pdat$yrse=yearfit$se

mscale=c(min(pdat$TPlwr, na.rm=TRUE)-sscal,
         max(pdat$TPupr, na.rm=T)+sscal)
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:3){
  
   datre=pdat[pdat$Station_Acronym==stvec[i],]
  plot(metric~year, xaxt="n",pch=16, data=datre,
       ann=FALSE, xlim=xyear,ylim=mscale, las=2)
  lines(TPfit~year, lwd=2,col=i, data=datre)
  lines(TPupr~year, data=datre,col="purple")
  lines(TPlwr~year, data=datre,col="purple")
  x=c(rev(datre$TPlwr), datre$TPupr)
  y=c(rev(datre$year), datre$year)
  polygon(y,x, col=adjustcolor("grey",alpha.f=0.1), border=NA )
}


par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
derivatives(fm, interval="simultaneous")
dtall=as.data.frame(derivatives(fm, interval="simultaneous"))
dtall=as.data.frame(derivatives(fm, interval="simultaneous"))
dtall$smooth=as.factor(dtall$smooth)
dlit=split(dtall, dtall$smooth)
for (i in 2:4) {
dat=dlit[i]
dat=as.data.frame(dlit[i])
colnames(dat)=colnames(dtall)
plot(derivative~data, data=dat, ylim=c(-1,1), typ="l", lwd=2)
lines(lower~data, data=dat)
lines(upper~data, data=dat)
abline(h=0)
abline(v=1994)

}
