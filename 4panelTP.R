# 4 panel of TP effects either full or partial smooths

#fm is light model with partial (fmvec[7])
#fc is chlorphy model with partial (fcvec[7])
# fmvec[12] light model full
# fcvec[12] chl model full



#plot TP smooths for light


#####

#tiff('baseR_figure.tiff', width =130, 
#     height = 80, pointsize = 12,
#     units = 'mm', res = 300)

jpeg('baseR_figure.jpeg', width =150,height = 100,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

par(mfcol=c(2,2), mar=c(0.5,6,1,1),oma=c(4,1,0,0),cex=.8,
#, oma=c(4,3,0,0),
    cex.axis=.9, cex.lab=1.1, cex.main=1.2, cex.sub=1)
#ylim2=c(-.5,.55)


fmTP<- gam(as.formula(frmvec[12]) ,
           data = pdat, weights=weights, 
           
           na.action=na.omit,
           niterPQL=100,
           #correlation=corAR1(form = ~ year|Station_Acronym),
           
           method = "REML")
TPsoutfull=smooth_estimates(fmTP)


plot(pdat[,"metric"]~pdat$TP, pch=as.numeric(pdat$Station_Acronym)+14,
     col=pdat$Station_Acronym, xaxt="n", xlab="",
    
     data=pdat,bty='l', cex=.8, las=1, 
     #ylim=c(0, 2.25),
    # xlab=expression(Total~phosophorus~(mg~L^-1)),
     ylab=expression(Light~attenuation~(m^-1))) 

TPsoutfull$pred=TPsoutfull$est+summary(fmTP)$p.table[1]
lines(pred~TP, data=TPsoutfull,
      lwd=2,col="dark blue")

TPsoutfull$upper =  TPsoutfull$pred + 2*TPsoutfull$se
TPsoutfull$lower =  TPsoutfull$pred - 2*TPsoutfull$se
lines(upper~TP, data=TPsoutfull,col="dark blue")
lines(lower~TP, data=TPsoutfull,col="dark blue")
x=c(rev(TPsoutfull$lower), TPsoutfull$upper)
y=c(rev(TPsoutfull$TP), TPsoutfull$TP)
polygon(y,x, col=adjustcolor("dark blue",alpha.f=0.05), border=NA )

legend("topleft","(a)", bty="n")
  #     inset=c(-0.025,-0.07))

legend("bottomright",c("Belleville", "Hay Bay", "Conway"),
       bty="n",
       col=c(1,3,2), pch=c(15,17,16), cex=.9)
#     inset=c(-0.025,-0.07))




TPpout=partial_residuals(fm, select = "s(TP)", partial_match=TRUE)
TPsout=smooth_estimates(fm, smooth = "s(TP)", partial_match=TRUE)
pIdat=cbind(pdat,TPpout)



plot(pIdat[,"s(TP)"]~pIdat$TP, pch=as.numeric(pdat$Station_Acronym)+14,
     col=pdat$Station_Acronym, 
     #ylim=c(-0.3, 0.75),
     data=pIdat,bty='l', cex=.8, las=1, 
     xlab=expression(Total~phosophorus~(mg~L^-1)),
     ylab="Partial effect on \nlight attenutation")
lines(TPsout$est~TPsout$TP, lwd=2,col="dark blue")
# points(mgres~year, data=datsres, col="red")

TPsout$upper =  TPsout$est + 2*TPsout$se
TPsout$lower =  TPsout$est - 2*TPsout$se
lines(upper~TP, data=TPsout,col="dark blue")
lines(lower~TP, data=TPsout,col="dark blue")
x=c(rev(TPsout$lower), TPsout$upper)
y=c(rev(TPsout$TP), TPsout$TP)
polygon(y,x, col=adjustcolor("dark blue",alpha.f=0.05), border=NA )

legend("topleft","(b)", bty="n")
mtext(expression(Total~phosophorus~(mg~L^-1)), 
      at=0.33,
      side = 1, line = +2, outer = TRUE, 
      adj = NA, padj = NA, cex = .9, col = NA, font = NA)


##CHL

fcTP<- gam(as.formula(fcvec[12]) ,
           data = pcdat, weights=weights, 
           
           na.action=na.omit,
           niterPQL=100,
           #correlation=corAR1(form = ~ year|Station_Acronym),
           
           method = "REML")
cTPsoutfull=smooth_estimates(fcTP)


plot(pcdat[,"metric"]~pcdat$TP, 
     pch=as.numeric(pdat$Station_Acronym)+14,
     col=pdat$Station_Acronym, xaxt="n", xlab="",
     
     data=pcdat,bty='l', cex=.8, las=1, 
     ylim=c(0, 43),
      ylab=expression(Chlorophyll~a~(mu~g~L^-1))) 

pred=cTPsoutfull$est+summary(fcTP)$p.table[1]
lines(pred~cTPsoutfull$TP, lwd=2,col="dark blue")
# points(mgres~year, data=datsres, col="red")

cTPsoutfull$upper =  pred + 2*cTPsoutfull$se
cTPsoutfull$lower =  pred - 2*cTPsoutfull$se
lines(upper~TP, data=cTPsoutfull,col="dark blue")
lines(lower~TP, data=cTPsoutfull,col="dark blue")
x=c(rev(cTPsoutfull$lower), cTPsoutfull$upper)
y=c(rev(cTPsoutfull$TP), cTPsoutfull$TP)
polygon(y,x, col=adjustcolor("dark blue",alpha.f=0.05), border=NA )

legend("topleft","(c)", bty="n")
#     inset=c(-0.025,-0.07))





cTPpout=partial_residuals(fc, select = "s(TP)", partial_match=TRUE)
cTPsout=smooth_estimates(fc, smooth = "s(TP)", partial_match=TRUE)
pcTdat=cbind(pcdat,cTPpout)



plot(pcTdat[,"s(TP)"]~pcTdat$TP, pch=as.numeric(pcTdat$Station_Acronym)+14,
     col=pcTdat$Station_Acronym, 
     #ylim=c(-0.3, 0.75),
     data=pcTdat,bty='l', cex=.8, las=1, 
     xlab=expression(Total~phosophorus~(mg~L^-1)),
     ylab="Partial effect on \nchlorophyll a")
lines(cTPsout$est~cTPsout$TP, lwd=2,col="dark blue")
# points(mgres~year, data=datsres, col="red")

cTPsout$upper =  cTPsout$est + 2*cTPsout$se
cTPsout$lower =  cTPsout$est - 2*cTPsout$se
lines(upper~TP, data=cTPsout,col="dark blue")
lines(lower~TP, data=cTPsout,col="dark blue")
x=c(rev(cTPsout$lower), cTPsout$upper)
y=c(rev(cTPsout$TP), cTPsout$TP)
polygon(y,x, col=adjustcolor("dark blue",alpha.f=0.05),
        border=NA )

legend("topleft","(d)", bty="n")
mtext(expression(Total~phosophorus~(mg~L^-1)), 
      at=0.83,
      side = 1, line = +2, outer = TRUE, 
      adj = NA, padj = NA, cex = .9, col = NA, font = NA)

dev.off()