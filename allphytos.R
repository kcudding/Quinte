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
stvec=c("B", "HB",  "C")

stlabel=c("Belleville", "Hay Bay", "Conway")
#select sites on which to run analysis

par(mfrow=c(3,2), mar=c(0,4,0,3), oma=c(5,5,2,2),cex.axis=1.25)


spydat=agg[agg$Station_Acronym%in%stvec,]

xyear=c(1972,2008)
spydat=spydat[spydat$year%in%c(1972:2008),]
spydat$Station_Acronym=as.factor(spydat$Station_Acronym)
#str(spydat)






site=stvec
bxlab=c("","", "Breaks")
bylab=c("","BIC", "")
ly=expression(paste("Light Attentuated (", m^-1, ")"))
lylab=c("", ly, "")


#TP model
tp1=(lm(metric~TP+Station_Acronym, weights=weights, data=spydat))
tp2=(lm(metric~TP:Station_Acronym, weights=weights, data=spydat))
tp3=(lm(metric~TP:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
cop=AIC(tp1,tp2,tp3)
print(cop)
print(c(summary(tp1)$adj.r.squared,summary(tp2)$adj.r.squared,summary(tp3)$adj.r.squared))


#metric model

tp4=(lm(metric~TP+year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
tp5=(lm(metric~TP+year+Station_Acronym, weights=weights, data=spydat))
tp6=(lm(metric~TP:year:Station_Acronym+Station_Acronym, weights=weights, data=spydat))
tp7=(lm(metric~TP:Station_Acronym+year+Station_Acronym, weights=weights, data=spydat))
tp8=(lm(metric~TP+year:Station_Acronym, weights=weights, data=spydat))
tp9=(lm(metric~TP+year+Station_Acronym, weights=weights, data=spydat))
tp10=(lm(metric~TP:year:Station_Acronym, weights=weights, data=spydat))
tp11=(lm(metric~TP:Station_Acronym+year, weights=weights, data=spydat))





cop2=AIC(tp4,tp5,tp6, tp7, tp8, tp9, tp10, tp11)
print(cop2)
print(c(summary(tp4)$adj.r.squared,summary(tp5)$adj.r.squared,summary(tp7)$adj.r.squared,
        summary(tp8)$adj.r.squared))

m1=which.min(cop$AIC)
mod=row.names(cop)[m1]
bestmod=as.name(mod)

#summary(eval(bestmod))
mlab="phyto"

mpred=as.data.frame(predict(tp3, spydat, interval="confidence"))

spydat=cbind(spydat, mpred)
sscal=1
#spydat$tres=spydat$TP-spydat$tmodel
par(mfrow=c(3,1), mar=c(1,4,0,3), oma=c(5,5,2,2),
    cex.axis=1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:length(stvec)) {
  pdat=spydat[spydat$Station_Acronym==stvec[i],]
  
  mscale=c(min(pdat$metric, na.rm=TRUE)-sscal,
           max(pdat$metric, na.rm=T)+sscal)
 # mscale=c(0,4500)
  amscale=mscale
  mss=which(!complete.cases(pdat))
  if (length(mss)>0) pdat=pdat[-mss,]
  if (i==3) {
    
    amscale[2]=mscale[2]/2
  }
  ymax=round(max(pdat$metric, na.rm=T),2)
  plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
       ann=FALSE, xlim=xyear, ylim=amscale, las=2, yaxp=c(0,ymax,3))
  
  if ( i==3) {
     axis(side=1, las=1)
    mtext(mmetric,side=1)
  }
  
  legend("bottomleft",stlabel[i], bty="n",
         cex=1.5)
  
  lines(pdat$fit~pdat$year, col = 2, lwd=2)
  lines(pdat$lwr~pdat$year, col = "pink", lwd=2, lty=2)
  lines(pdat$upr~pdat$year, col = "pink", lwd=2, lty=2)
  abline(v=1995)
  abline(v=1988)
}

treg=summary(tp3)$coefficients

} 
#regco=kable(treg, digits=2,
#format="html", table.attr = "style='width:60%;'",
#   caption="TP~year")
#tsave=footnote(kable_classic(regco),
#         general = paste("adj. R squared=", round((summary(tp3)$adj.r.squared),2)),
#         footnote_as_chunk = TRUE,
#         general_title = "")
#tsave
#fname=paste0(mmetric,"Site", stvec[i],".png")
#print(tsave)
#cat('\n\n<!-- -->\n\n')



```

```{r, results='asis'}
reglist=list()
#tiff('baseR_figure.tiff', width =85, 
#     height = 180, pointsize = 12,
#     units = 'mm', res = 300)


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
  cat('\n\n<!-- -->\n\n')
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
    fm.riT2 <- lm(metric ~ time:year+TP+time, data=pdat, weights=weights)
    fm.riT <- lm(metric ~ time:year+TP, data=pdat, weights=weights)
    fm.ri <- lm(metric ~ time+TP+year, data=pdat, weights=weights)
    fm.riN <- lm(metric ~ TP+year, data=pdat,weights=weights)
    pAIC=(AIC(fm.ri, fm.riN,fm.riT, fm.riT2))
    
    print(pAIC)
    cat('\n\n<!-- -->\n\n')
    summary(fm.ri)
    fm.ri=fm.riT2
  } else {
    fm.ri <- lm(metric ~ TP+year, data=pdat,weights=weights)
    fm.riTP <- lm(metric ~ TP, data=pdat,weights=weights)
    fm.riyr <- lm(metric ~ year, data=pdat,weights=weights)
  }
  
  #table of output
  
  reglist[[i]]=fm.ri  
  
  summary(fm.ri)$adj.r.squared
  
  
  amscale=mscale
  amscale=c(0, amscale[2])
  if (i==3) amscale[2]=mscale[2]/2
  plot(metric~year, data=pdat, type="b",xaxt="n",pch=16,col="grey",
       ann=FALSE, xlim=xyear, ylim=amscale,las=2,yaxp=c(0,amscale[2],4))
  
  if ( i==3)axis(side=1, las=1)
  legend("topright",stlabel[i], bty="n",
         cex=1.5)
  
  bfm=predict( fm.ri, newdata=pdat)
  
  lines(bfm~pdat$year, col = 2, lwd=2)
  
  
  #convert confidence intervals from julian day to dates
  #and draw
  
  #draw breakpoints & confidence intervals 
  if (sum(!is.na(b$breakpoints))>0){
    cis=confint(b)$confint
    cll=cis[summary(fm.ri)$coefficients[2:(1+length(b$breakpoints)),4]<1.,]
    ci=cis[summary(fm.ri)$coefficients[2:(1+length(b$breakpoints)),4]<0.05,]
    if (length(nrow(cll))!= 0) {
      if (nrow(cll)== 1) abline(v=pdat$year[cll[2]], lty=2, col="light grey")
      if (nrow(cll)>1) abline(v=pdat$year[cll[,2]], lty=2, col="light grey")
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
      xn=pdat$year[ci[1]]
      xx=pdat$year[ci[3]]
      yx=yn
      segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
      text(pdat$year[ci[2]]-2,amscale[2]-.5, pdat$year[ci[2]])
      abline(v=pdat$year[ci[2]], lty=2, col="blue")
    } else {
      yn=rep(hgt,tm)
      
      s=1:tm
      
      xn=pdat$year[ci[,1]]
      xx=pdat$year[ci[,3]]
      yx=yn
      segments(xn[s], yn[s], xx[s], yx[s], col="blue", lwd=2, lend="square")
      abline(v=pdat$year[ci[,2]], lty=2, col="blue")
      text(pdat$year[ci[,2]]-2,amscale[2]-.5, pdat$year[ci[,2]])
    }
  }
  if (i==1)legend("bottom", c("data", "model", "break",
                              "break CI", "ns break"), 
                  lty=c(1,1,2,1,2), col=c("grey", "red","blue", "blue", "grey"), ncol=2,xjust=1, yjust=1,y.intersp = .75,
                  pch=c(16,NA,NA,NA, NA),
                  inset=c(0, -.025),
                  bty="n", cex=1.2)
}

mtext(mlab, side = 2, line = +1, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)
mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)



dev.off()

##########END BREAKPOINT
```


```{r}
i=1
regs=as.data.frame(summary(reglist[[i]])$coefficients) 
sigstar=ifelse (regs[,4]<0.05, "*", "") 
regs[,4]=formatC(regs[,4], format = "e", digits = 1)
regs[,4]=paste(regs[,4],sigstar)

rsq=round((summary(reglist[[i]])$adj.r.squared),4)


regco=knitr::kable((regs), digits=4,
                   #format="html", table.attr = "style='width:60%;'",
                   caption=paste(mmetric,"Site", stvec[i]))
#tsave=footnote(kable_classic(regco),
#        general = paste("adj. R squared=", #round((summary(reglist[[i]])$adj.r.squared),4)),
#        footnote_as_chunk = TRUE,
#       general_title = "")
regco
print(paste("adj.r.squared", rsq))
#fname=paste0(mmetric,"Site", stvec[i],".png")
#tsave
#cat('\n\n<!-- -->\n\n')

#save_kable(tsave,file=fname, density=600, zoom = 1.25, bs_theme = "flatly")


```

```{r}
i=2
regs=as.data.frame(summary(reglist[[i]])$coefficients) 
sigstar=ifelse (regs[,4]<0.05, "*", "") 
regs[,4]=formatC(regs[,4], format = "e", digits = 1)
regs[,4]=paste(regs[,4],sigstar)
rsq=round((summary(reglist[[i]])$adj.r.squared),4)
regco=knitr::kable((regs), digits=4,
                   #format="html", table.attr = "style='width:60%;'",
                   caption=paste(mmetric,"Site", stvec[i]))
#tsave=footnote(kable_classic(regco),
#        general = paste("adj. R squared=", #round((summary(reglist[[i]])$adj.r.squared),4)),
#        footnote_as_chunk = TRUE,
#       general_title = "")

#fname=paste0(mmetric,"Site", stvec[i],".png")
#tsave
#cat('\n\n<!-- -->\n\n')
regco
print(paste("adj.r.squared", rsq))
#save_kable(tsave,file=fname, density=600, zoom = 1.25, bs_theme = "flatly")


```

```{r}
i=3
regs=as.data.frame(summary(reglist[[i]])$coefficients) 
sigstar=ifelse (regs[,4]<0.05, "*", "") 
regs[,4]=formatC(regs[,4], format = "e", digits = 1)
regs[,4]=paste(regs[,4],sigstar)
rsq=round((summary(reglist[[i]])$adj.r.squared),4)
regco=knitr::kable((regs), digits=4,
                   #format="html", table.attr = "style='width:60%;'",
                   caption=paste(mmetric,"Site", stvec[i]))
#tsave=footnote(kable_classic(regco),
#        general = paste("adj. R squared=", #round((summary(reglist[[i]])$adj.r.squared),4)),
#        footnote_as_chunk = TRUE,
#       general_title = "")

#fname=paste0(mmetric,"Site", stvec[i],".png")
#tsave
#cat('\n\n<!-- -->\n\n')
regco
print(paste("adj.r.squared", rsq))
#save_kable(tsave,file=fname, density=600, zoom = 1.25, bs_theme = "flatly")


```
