temp0=read.csv("Trentontemps_1964_04_19_to_1991_09_04.csv")

temp1=read.csv("Trentontemps_1991_09_05_to_2019_03_26.csv")
temps=rbind(temp0,temp1)
date_t=as.POSIXct(temps$LOCAL_DATE)
mean_t=temps$MEAN_TEMPERATURE

plot(mean_t~date_t, data=temps, type="l")

tempsmonth=temps[temps$LOCAL_MONTH%in%c(3,4,5),]
temps=tempsmonth
#ann_mean=aggregate(temps$MEAN_TEMPERATURE~temps$LOCAL_YEAR+LOCAL_MONTH, data=temps, FUN = "mean" )

#ann_mean$date=as.Date(paste0("01", "/", ann_mean$LOCAL_MONTH, "/", ann_mean$`temps$LOCAL_YEAR`), format="%d/%m/%Y")

ann_mean=aggregate(temps$MEAN_TEMPERATURE~temps$LOCAL_YEAR, data=temps, FUN = "mean" )
ann_min=aggregate(temps$MIN_TEMPERATURE~temps$LOCAL_YEAR, data=temps, FUN = "mean" )
ann_max=aggregate(temps$MAX_TEMPERATURE~temps$LOCAL_YEAR, data=temps, FUN = "mean" )
ann_wgt=aggregate(temps$MEAN_TEMPERATURE~temps$LOCAL_YEAR, data=temps,FUN=length)
colnames(ann_wgt)[2]="n_samp"
anns=merge(ann_mean, ann_min)
anns=merge(anns,ann_max)
anns=merge(anns, ann_wgt)

colnames(anns)=c("year", "mean", "avg_min", "avg_max", "n")

library("strucchange")

temp_comm=merge(spydat, anns)
summary(lm(epar.mean~mean, weights=epar.wgt, data=temp_comm))





spydat=agg[agg$Station_Acronym%in%stvec,]
xyear=c(1979,2015)
an_sel=anns[anns$year%in%c(1986:2007),]
spydat=spydat[spydat$year%in%c(1986:2007),]
spydat$Station_Acronym=as.factor(spydat$Station_Acronym)



plot(mean~year, data=an_sel, type="b")
lines(avg_min~year, data=an_sel, type="l", lty=2)
lines(avg_max~year, data=an_sel, type="l", lty=2)


tp1=lm(spydat$Chrysophyceae.mean~spydat$TP)
tp2=lm(spydat$Chrysophyceae.mean~spydat$TP+an_sel$mean)
     AIC(tp1,tp2)
     summary(tp1)



spyyear=spydat[complete.cases(spydat),]



