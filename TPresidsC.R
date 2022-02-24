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
xyear=c(1972, 2008)
spydat=agg[agg$Station_Acronym%in%stvec,]

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)
pdat=spydat
pdat=pdat[pdat$year%in%c(1972:2008),]
xyear=c(1972,2008)
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
frcvec1=c(form1,form2, form3, form4)
sc=vector()
for (fr in 1:length(frcvec1)){
  mc<- gam(as.formula(frcvec1[fr]) ,
          data = pdat, weights=weights, 
          na.action=na.omit,
          niterPQL=100,
          #correlation=corAR1(form = ~ year|Station_Acronym),
          method = "REML")
  sc[fr]=AIC(mc)
  
}

ind=which.min(sc)
ind=3
mc<- gam(as.formula(frcvec1[ind]) ,
        data = pdat, weights=weights, 
        na.action=na.omit,
        niterPQL=100,
        #correlation=corAR1(form = ~ year|Station_Acronym),
        method = "REML")
summary(mc)
AIC(mc)
concurvity(mc, full=TRUE)
round(concurvity(mc, full=FALSE)$worst,2)

draw(mc, residuals=TRUE, rug=FALSE)
TPfit=predict(mc, newdata=pdat, se.fit=TRUE)
pdat$TPfit=TPfit$fit
pdat$TPfitse=TPfit$se
pdat$TPres=pdat$metric-pdat$TPfit
crit.t <- qt(0.975, df = df.residual(mc)) 
pdat$TPupr = TPfit$fit + (crit.t *TPfit$se)
pdat$TPlwr = TPfit$fit + (crit.t *TPfit$se)

mscale=c(min(newYear$lower, na.rm=TRUE)-sscal,
         max(pdat$metric, na.rm=T)+sscal)




form6="TPres ~  Station_Acronym+s(year,by=Station_Acronym, bs='tp', k=20,m=1)" 
form7="TPres ~  s(year,by=Station_Acronym,bs='tp',k=20,m=2)"
form8="TPres ~  s(year,bs='tp',k=20,m=2)"




frcvec=c(form6,form7, form8)
rm(mm)
mm=list()
sc=vector()
for (fr in 1:length(frcvec)){
 frc<- gam(as.formula(frcvec[fr]) ,
          data = pdat, weights=weights, 
          
          na.action=na.omit,
          niterPQL=100,
          #correlation=corAR1(form = ~ year|Station_Acronym),
          method = "REML")
  sc[fr]=AIC(frc)
  
  
  
}
ind=which.min(sc)
ind=1
frc<- gam(as.formula(frcvec[ind]) ,
         data = pdat, weights=weights, 
         
         na.action=na.omit,
         niterPQL=100,
         correlation=corAR1(form = ~ year|Station_Acronym),
         
         method = "REML")

summary(frc)
AIC(frc)
concurvity(frc, full=TRUE)
round(concurvity(frc, full=FALSE)$worst,2)
draw(m, residuals=TRUE)
draw(frc, residuals=TRUE, rug=FALSE)

plot(frc, seWithMean=TRUE, shift=coef(frc)[1], rug=FALSE, 
     shade=TRUE,residuals=TRUE, pch=16)
#by.resids=TRUE)

frcfit=fitted_values(frc)
pdat=cbind(pdat,frcfit)

panlab=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(j)", "(i)" )
jpeg('baseR_figure.jpeg', width =130,height = 80,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

par(mfcol=c(3,2), mar=c(0.5,4,1,1), oma=c(4,1,0,0),
    cex.axis=1.1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:3){
  ylim1=c(-15,15)
 # if (i==3) ylim1=ylim1/2
  datre=pdat[pdat$Station_Acronym==stvec[i],]
  plot(TPres~year, xaxt="n",pch=16, data=datre,bty='l',cex=.8,
       ann=FALSE,  xlim=c(1972.5,2007.5),
       ylim=ylim1,
       las=2, col="dark green")
  lines(fitted~year, lwd=2,data=datre,col="dark green")
  # lines(upper~year, data=datre,col="purple")
  #  lines(lower~year, data=datre,col="purple")
  x=c(rev(datre$lower), datre$upper)
  y=c(rev(datre$year), datre$year)
  polygon(y,x, col=adjustcolor("dark green",alpha.f=0.1), border=NA )
  legend("topright",panlab[i], bty="n",inset=c(0, -.15),
         cex=1.)
  #  text(y=35, x=1990,stlabel[i], bty="n",
  #      cex=1.25 )
  # title(stlabel[i], cex=0.8)
  mtext(stlabel[i], side = 3, line = 0, outer = FALSE, at = NA,
        adj = NA, padj = NA, cex = 0.8, col = NA, font = NA)
  if ( i==3)axis(side=1, las=1)
}
mtext("residuals of chlorophyll a", side = 2, line = -1.0, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)


dtall=as.data.frame(derivatives(frc,n=37, interval="simultaneous"))

dtall$smooth=as.factor(dtall$smooth)
dlit=split(dtall, dtall$smooth)
yaxl=c(-5.5,5.5)
dseq=c(1,3,2)

for (j in seq_along(dseq)) {
  if (j==3) yaxl=yaxl/10
  i=dseq[j]
  df.new = dlit[[i]][seq(1, nrow(dlit[[i]]), 5), ]
  dat=as.data.frame(df.new)
  colnames(dat)=colnames(dtall)
  dat=dtall[dtall$smooth==levels(dtall$smooth)[i],]
  
  plot(derivative~data, data=dat,xaxt="n",ann=FALSE,bty='l',las=2,
       ylim=yaxl, typ="l", lwd=2, col="dark red", xlim=c(1972.5,2007.5))
  lines(lower~data, data=dat, col="dark red")
  lines(upper~data, data=dat,col="dark red")
  x=c(rev(dat$lower), dat$upper)
  x[x>yaxl[2]]<-yaxl[2]+.1
  x[x< yaxl[1]]<-yaxl[1]-.1
  y=c(rev(dat$data), dat$data)
  polygon(y,x, col=adjustcolor("dark red",alpha.f=0.1), border=NA )
  
  abline(h=0, col="dark grey", lwd=2)
  abline(v=1994, col="dark grey", lwd=2)
  if (j==3) text(1990, 0.4, "1994", cex=0.9)
  if (j==3)axis(side=1, las=1)
  if (j==2) mtext("Derivative of smooth", side = 2, line = 3., outer = FALSE, at = NA,
                  adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
  legend("topright",panlab[j+3], bty="n",
         cex=1., inset=c(0, -.15))
  
}
dev.off()
