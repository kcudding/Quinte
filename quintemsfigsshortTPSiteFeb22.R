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
xyear=c(1972,2015)
spydat=agg[agg$Station_Acronym%in%stvec,]

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)
pdat=spydat
pdat=pdat[pdat$year%in%c(1972:2015),]
xyear=c(1972,2015)
str(pdat)
N=1000 #number of points at which to evaluate the smooth

pdat$Station_Acronym=as.factor(pdat$Station_Acronym)

form1="metric ~  s(TP, bs='tp')+s(Station_Acronym, bs='re')"  #global
form2="metric ~  s(TP, bs='tp')+s(TP, Station_Acronym,bs='fs',m =2)" #gs
form3="metric ~  s(TP, bs='tp')+s(TP,by=Station_Acronym, 
bs='tp', m=1)+s(Station_Acronym, bs='re')" #gi
form4="metric ~  s(TP, Station_Acronym,bs='fs',m =2)" #s
form5="metric ~  s(TP,by=Station_Acronym, 
bs='tp', m=1)+s(Station_Acronym, bs='re')" #i
form6="metric ~  s(TP, bs='tp')"
frmvec=c(form1, form2,form3, form4, form5, form6)

sc=vector()
for (fr in 1:length(frmvec)){
m<- gam(as.formula(frmvec[fr]) ,
          data = pdat, weights=weights, 
          
          na.action=na.omit,
          niterPQL=100,
#correlation=corAR1(form = ~ year|Station_Acronym),
method = "REML")
sc[fr]=AIC(m)

}

ind=which.min(sc)
fm<- gam(as.formula(frmvec[ind]) ,
        data = pdat, weights=weights, 
        
        na.action=na.omit,
        niterPQL=100,
        #correlation=corAR1(form = ~ year|Station_Acronym),
        method = "REML")





ad=as.data.frame(derivatives(fm, term="s(year):Station_AcronymB"))

par(mfrow=c(3,1))
vis.gam(fm, view=c('year','TP'), cond=list(Station_Acronym="C"), 
        type='response', theta=30, main="Station C", zlab=mmetric)
vis.gam(fm, view=c('year','TP'), cond=list(Station_Acronym="HB"), 
        type='response', theta=30, main="Station HB", zlab=mmetric)
vis.gam(fm, view=c('year','TP'), cond=list(Station_Acronym="B"), 
        type='response', theta=30, main="Station B", zlab=mmetric)

## create new data to predict at; 200 evenly-spaced values over `Year` 

newYear<-with(pdat,expand.grid(year = seq(from=min(year, na.rm=TRUE),to=max(year, na.rm=TRUE),  
            length.out = 200), 
            TP = seq(from=min(TP, na.rm=TRUE), to=max(TP, na.rm=TRUE), length.out=200),
            Station_Acronym = c("B","HB", "C")))


## Predict from the fitted model; note we predict from the $gam part
newYear <- cbind(newYear, data.frame(predict(fm, newYear, se.fit = TRUE)))
newYear <- cbind(pdat, data.frame(predict(fm, pdat, se.fit = TRUE)))
## Create the confidence interval 
crit.t <- qt(0.975, df = df.residual(fm)) 
newYear <- transform(newYear,
upper = fit + (crit.t * se.fit), lower = fit - (crit.t * se.fit))

mscale=c(min(newYear$lower)-sscal,
         max(pdat$metric, na.rm=T)+sscal)
ad=as.data.frame(derivatives(fm, term="s(year):Station_AcronymB"), interval="simultaneous")
plot(x=ad[,3], y=ad[,4])
lines(x=ad$data, y=ad$upper, lty=2)
lines(x=ad$data, y=ad$lower, lty=2)
abline(h=0)
abline(v=1994)

for (i in 1:3){
  
  datpl=newYear[newYear$Station_Acronym==stvec[i],]
  datre=pdat[pdat$Station_Acronym==stvec[i],]
  plot(metric~year, xaxt="n",pch=16, data=datre,
       ann=FALSE, xlim=xyear,ylim=mscale, las=2)
  lines(datpl$fit~datpl$year, lwd=2,col=i)
  lines(upper~year, data=datpl,col="purple")
  lines(lower~year, data=datpl,col="purple")
  x=c(rev(datpl$lower), datpl$upper)
  y=c(rev(datpl$year), datpl$year)
  polygon(y,x, col=adjustcolor("grey",alpha.f=0.1), border=NA )
}

nsim=50
small.d <- fderiv(fm, newdata = newYear, n = N) 
small.sint <- with(newYear,
cbind(confint(small.d, nsim = nsim, type = "simultaneous"),
                         Year = year))

rscale=c(min(small.sint$lower)-sscal,sscal)
ticks <- round(seq(from=rscale[1], to=rscale[2], length=6),2)

  if (mmetric=="epar" & i==2) {
    ticks=c(-0.05, -0.03, -0.01, 0.00,  0.01)
    rscale=c(min(small.sint$lower),sscal)
  }



plot(small.sint$est~small.sint$Year,ann=FALSE,xaxt = "n",
     type="l", lwd=2,xaxt="n",col="dark green",yaxt="n",
     ylim=rscale,xlim=xyear, las=2)


if ( i==3) {axis(side=1, las=1)}
if (i==2) {mtext(rlab,cex=1.3,line=4,side=2 ) }
  
axis(side=2,  las=1, at=ticks, labels=ticks )

lines(lower~Year, data=small.sint, col="green")
lines(upper~Year, data=small.sint, col="green")

x=c(rev(small.sint$lower), small.sint$upper)
y=c(rev(small.sint$Year), small.sint$Year)
polygon(y,x, col=adjustcolor("green",alpha.f=0.1), border=NA )

abline(h=0, lwd=2, col=adjustcolor("dark gray",alpha.f=0.75))
abline(v=1994, col="red", lwd=2)
#legend("bottomright",stlabel[i], bty="n",cex=1.5)
}


#mtext(expression(Light~attentuation~(m^-1)),side=2,line=0,outer=TRUE,cex=1.3,las=0)
mtext(mlab,side=2,line=0,outer=TRUE,cex=1.3,las=0)

mtext("Year",side=1,line=3,outer=TRUE,cex=1.3)


