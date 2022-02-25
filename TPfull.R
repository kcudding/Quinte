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

#mlab=expression(Light~attenuation~(m^-1)) 
mlab=expression(Chlorophyll~a~(mu~g~L^-1))
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
xyear=c(1972,2008)
spydat=agg[agg$Station_Acronym%in%stvec,]

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)
pdat=spydat
pdat=pdat[pdat$year%in%c(1972:2008),]
xyear=c(1972,2008)
str(pdat)
N=1000 #number of points at which to evaluate the smooth

pdat$Station_Acronym=as.factor(pdat$Station_Acronym)
pcdat=pdat

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)

str(pdat)
N=1000 #number of points at which to evaluate the smooth


form1="metric ~  te(TP,year,  bs='tp')"
form2="metric ~  te(TP,year, bs='tp', m=2)+t2(TP,year,Station_Acronym, bs=c('tp','tp','re'),  m=2, full=TRUE)"  #gs
form3="metric ~  Station_Acronym+te(TP,year, bs='tp', m=2)+te(TP,year,by=Station_Acronym, bs='tp', m=1)" #gi
form4="metric ~  t2(TP,year,Station_Acronym, bs=c('tp','tp','re'),m=2, full=TRUE)" #s
form5="metric ~  Station_Acronym+te(TP,year,by=Station_Acronym, bs=c('tp', 'tp'), m=2)" #i
form6="metric ~  Station_Acronym+s(TP,bs='tp',m=2)+te(TP,year,by=Station_Acronym,  bs='tp',  m=1)" 
form7="metric ~  Station_Acronym+s(TP,bs='tp',m=2)+s(year,by=Station_Acronym, bs='tp', m=2)" 
form8="metric ~  Station_Acronym+s(TP,bs='tp',m=2)+s(year,bs='tp',m=2)" 
form9="metric ~  s(TP,bs='tp',m=2)+s(year,by=Station_Acronym,bs='tp',m=2)"
form10="metric ~  Station_Acronym+TP+s(year,by=Station_Acronym,bs='tp',m=2)"
form11="metric ~  TP+s(year,by=Station_Acronym,bs='tp',m=2)"
form12="metric ~  s(TP)"
fcvec=c(form1,form2,form3, form4, form5, form6, form7, form8, 
         form9, form10, form11, form12)
rm(mm)
mm=list()
sc=vector()

for (fr in 1:length(fcvec)){
  print(fr)
  m<- gam(as.formula(fcvec[fr]) ,
          data = pdat, weights=weights, 
          
          na.action=na.omit,
          niterPQL=100,
          #correlation=corAR1(form = ~ year|Station_Acronym),
          method = "REML")
  sc[fr]=AIC(m)
  
  
  
}
sc
scc=sc
sc-min(sc)
ind=which.min(sc)
ind=7
fc<- gam(as.formula(fcvec[ind]) ,
         data = pdat, weights=weights, 
         
         na.action=na.omit,
         niterPQL=100,
         #correlation=corAR1(form = ~ year|Station_Acronym),
         
         method = "REML")

summary(fc)
AIC(fc)
concurvity(fc, full=TRUE)
round(concurvity(fc, full=FALSE)$worst,2)


redform7="metric ~  Station_Acronym+s(year,by=Station_Acronym, bs='tp', m=2)" 

fcred<- gam(as.formula(fcvec[12]) ,
            data = pdat, weights=weights, 
            
            na.action=na.omit,
            niterPQL=100,
            # correlation=corAR1(form = ~ year|Station_Acronym),
            
            method = "REML")

summary(fcred)




draw(fc, residuals=TRUE, rug=FALSE)

drt=plot(fc, seWithMean=TRUE, shift=coef(fc)[1], rug=FALSE, 
     shade=TRUE,residuals=TRUE, pch=16)
#by.resids=TRUE)
dort=draw(fc, rug=FALSE, residuals=TRUE)

dort+labs(title = 'Main title', 
          subtitle = 'My subtitle', caption = 'My caption')
draw(m, residuals=TRUE, rug=FALSE)
fcfit=fitted_values(fc)
pdat=cbind(pdat,fcfit)

fcderiv=as.data.frame(derivatives(fc, interval="simultaneous"))
fcderiv$smooth=as.factor(fcderiv$smooth)
lfc=split(fcderiv, fcderiv$smooth)
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,4,2,2),
    cex.axis=0.8, cex.lab=1.75, cex.main=1.2, cex.sub=1)
dseq=c(1,3,2)
for (i in seq_along(dseq)) {
  dat=as.data.frame(lfc[i])
  colnames(dat)=colnames(fcderiv)
  plot(derivative~data, data=dat, ylim=c(-1.1,1.1), 
       typ="l", lwd=2)
  lines(lower~data, data=dat)
  lines(upper~data, data=dat)
  abline(h=0)
  abline(v=1994)
}

smooth_estimates(fc)
yearfit=predict(fc, newdata=pdat, se.fit=TRUE)
pdat$yrfit=yearfit$fit
pdat$yrse=yearfit$se

mscale=c(0,
         max(pdat$upper, na.rm=T)+sscal)
#####

#tiff('baseR_figure.tiff', width =130, 
#     height = 80, pointsize = 12,
#     units = 'mm', res = 300)

jpeg('baseR_figure.jpeg', width =130,height = 80,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

panlab=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(j)", "(i)" )


par(mfcol=c(3,3), mar=c(0.5,4,1,1), oma=c(4,1,0,0),
    cex.axis=1.1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:3){
  ylim1=c(0,43)
  if (i==3) ylim1=ylim1/2
  datre=pdat[pdat$Station_Acronym==stvec[i],]
  plot(metric~year, xaxt="n",pch=16, data=datre,bty='l',cex=.8,
       ann=FALSE,  xlim=c(1972.5,2007.5),ylim=ylim1,las=2)
  lines(fitted~year, lwd=2,data=datre,col=1)
# lines(upper~year, data=datre,col="purple")
#  lines(lower~year, data=datre,col="purple")
  x=c(rev(datre$lower), datre$upper)
  y=c(rev(datre$year), datre$year)
  polygon(y,x, col=adjustcolor("grey",alpha.f=0.3), border=NA )
  legend("topright",panlab[i], bty="n",
         cex=1.)
#  text(y=35, x=1990,stlabel[i], bty="n",
 #      cex=1.25 )
 # title(stlabel[i], cex=0.8)
  mtext(stlabel[i], side = 3, line = 0, outer = FALSE, at = NA,
        adj = NA, padj = NA, cex = 0.8, col = NA, font = NA)
  if ( i==3)axis(side=1, las=1)
}
mtext(mlab, side = 2, line = -1.5, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)

pout=partial_residuals(fc, select = "s(year)", partial_match=TRUE)
sout=smooth_estimates(fc, smooth = "s(year)", partial_match=TRUE)
pout=pout[,c(1,3,2)]
# or selected smooths
dfsm=edf(fc, smooth = colnames(pout))
dfsm=as.data.frame(dfsm)
drt[[2]]$p.resid

colnames(pout)=c("s(year):B", "s(year):HB", "s(year):C")

pdat=cbind(pdat,pout)
pdat$mgB=drt[[2]]$p.resid
pdat$mgHB=drt[[2]]$p.resid
#par(mfcol=c(3,2), mar=c(0.5,3,1,1), oma=c(4,3,0,0),
 #   cex.axis=1.2, cex.lab=1.75, cex.main=1.2, cex.sub=1)

for (i in 1:3){
  mgres=drt[[i+1]]$p.resid
  wind=which(pdat$Station_Acronym==stvec[i])
  mgres=mgres[wind]
  datsm=sout[sout$Station_Acronym==stvec[i],]
  datsres=pdat[wind,]
  yname=(paste0("s(year):",stvec[i]))
 
  plot(datsres[,yname]~year, xaxt="n",pch=16,col="grey", 
       data=datsres,bty='l', cex=.8,
       ann=FALSE, xlim=c(1972.5,2007.5),ylim=c(-15,15),las=2)
  lines(est~year, lwd=2,data=datsm,col="dark green")
 # points(mgres~year, data=datsres, col="red")
  datsm$upper = datsm$est + (dfsm[i,2] * datsm$se)
  datsm$lower = datsm$est - (dfsm[i,2] * datsm$se)
  datsm$upper =  datsm$est + 2*datsm$se
  datsm$lower =  datsm$est - 2*datsm$se
   lines(upper~year, data=datsm,col="dark green")
    lines(lower~year, data=datsm,col="dark green")
  x=c(rev(datsm$lower), datsm$upper)
  y=c(rev(datsm$year), datsm$year)
  polygon(y,x, col=adjustcolor("dark green",alpha.f=0.1), border=NA )
  
  if (i==2) mtext("Partial effect", side = 2, line = 2.75, outer = FALSE, at = NA,
                adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
  if ( i==3)axis(side=1, las=1)
  #legend("topright",stlabel[i], bty="n",cex=1.5)
  legend("topright",panlab[i+3], bty="n",
         cex=1.)
}


derivatives(fc, interval="simultaneous")
dtall=as.data.frame(derivatives(fc, interval="simultaneous"))
dtall=as.data.frame(derivatives(fc, interval="simultaneous"))
dtall$smooth=as.factor(dtall$smooth)
dlit=split(dtall, dtall$smooth)
yaxl=c(-1.5,1.5)
dseq=c(2,4,3)

for (j in seq_along(dseq)) {
  i=dseq[j]
  dat=dlit[i]
  dat=as.data.frame(dlit[i])
  colnames(dat)=colnames(dtall)
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
  if (i==3) text(1988, 1.4, "1994", cex=0.9)
  if (i==3)axis(side=1, las=1)
  if (i==4) mtext("Derivative of smooth", side = 2, line = 3.25, outer = FALSE, at = NA,
                  adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
  legend("topright",panlab[i+6], bty="n",
         cex=1.)
  
}
dev.off()
