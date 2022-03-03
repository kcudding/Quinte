#using Sim.DiffProc and fitsde
rm(list=ls())
library("Sim.DiffProc")
#install.packages("strucchange")
library("strucchange")

#install.packages("mgcv")
library("mgcv")
#install.packages("devtools") 
#library("devtools")
#devtools::install_github("gavinsimpson/gratia")
library("gratia")

#parameters
a= 0.64254962
#a=.779
p=12
h=1.0
r=.6
b=1.0
sigma=0.05
start=1.13
#start=0.83

library(rootSolve)
# the model : lake

#create simulated data
f <- expression((r*x^p)/(x^p+h^p)+a-b*x)
g <- expression(sigma)

rate <- function(x, a = 0.67) (r*x^p)/(x^p+h^p)+a-b*x
# find all roots within the interval [0,10]
Eq   <- uniroot.all(rate, c(0, 2))
start=Eq[3]

ap=seq(0.64254962, 0.55, -0.01)
ap=0.6425
#ap=0.5025496
reps=1
freqT=40*12
lenT=40
tcon=lenT/freqT
dval=0.01
panlab=c("(a)","(b)", "(c)", "(d)", "(e)", "(f)")

par(mfrow=c(3,2), mar=c(0.5,7,1,1), oma=c(4,1,1,1),cex.axis=1.25)

#par(mfcol=c(3,3), mar=c(0.5,4,1,1), oma=c(4,1,0,0),
 #   cex.axis=1.1, cex.lab=1.75, cex.main=1.2, cex.sub=1)
for (i in 1:3) {
tr = matrix(nrow=reps,ncol=length(ap))
eqd=vector()
lg=vector()


for(j in 1:length(ap)){
 
  print(ap[j])
  a = ap[j]
rate <- function(x, a = ap[j]) (r*x^p)/(x^p+h^p)+a-b*x
EqT   <- uniroot.all(rate, c(0, 2))
ed=EqT[1]

eqd[j]=ed
print(ed)
#TEX.sde(object=c(drift = f, diffusion = g))
f <- expression((r*x^p)/(x^p+h^p)+a-b*x)
if (i==1) sigma=.01
if (i==2) sigma=.05
if (i==3) sigma=.1
  #g=expression(sigma)


for (k in 1:reps) {
HWV <- snssde1d(drift=f,diffusion=g,x0=start,N=freqT,T=lenT)
#plot(HWV$X, col=i)

tr[k,j]=min(which((abs(HWV$X-ed)<dval)))
}
}
samp=HWV$X
#sampd=data.frame(X=HWV$X)
#samp$time <- 1:nrow(samp)
sdat=ts(samp, start=1978, end=2018, freq=12)
ot<-aggregate.ts(sdat, nfrequency =1, mean)


yr=seq(from=1978, to=2017)
dotx=as.vector(ot)
dot=data.frame(metric=dotx, year=yr)
N=1000 #number of points at which to evaluat the smoot
m <- gamm(metric ~ s(year, k=40), data=dot, method = "REML")
## create new data to predict at; 200 evenly-spaced values over `Year` 
newYear <- with(dot, data.frame(year = seq(min(year), max(year), length.out = 200))) 
## Predict from the fitted model; note we predict from the $gam part
newYear <- cbind(newYear, data.frame(predict(m$gam, newYear, se.fit = TRUE)))
## Create the confidence interval 
crit.t <- qt(0.975, 
             df = df.residual(m$gam)) 
newYear <- transform(newYear,
                     upper = fit + (crit.t * se.fit), lower = fit - (crit.t * se.fit))
plot(metric~year, data=dot, type="b",
     bty='l',cex=.8, xaxt="n", las=2,
     ann=FALSE,ylim=c(0.5,1.3))
#plot(metric~year, data=pdat, type="b", main=paste("light extinction coef",stvec[i] ))
#abline(v=1994, col="red", lwd=2)
if (i==2) mtext("Bioturbidity", side = 2, line = -2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
if ( i==3)axis(side=1, las=1)


lines(newYear$fit~newYear$year, lwd=3)
lines(upper~year, data=newYear,col="green")
lines(lower~year, data=newYear,col="green")
legend("topright", panlab[i], bty="n", cex=1.2) 
       #inset=c(-.2,0))



nsim=50
small.d <- fderiv(m, newdata = newYear, n = N) 
small.sint <- with(newYear,
                   cbind(confint(small.d, nsim = nsim, type = "simultaneous"),
                         Year = year))



plot(small.sint$est~small.sint$Year,bty="l",ylim=c(-.1,.1),
     type="l", lwd=2, xaxt="n", las=2,ann=FALSE)
lines(lower~Year, data=small.sint, col="green")
lines(upper~Year, data=small.sint, col="green")
abline(h=0)
#abline(v=1994, col="red", lwd=2)
if (i==2) mtext("Derivative of smooth", side = 2, line = +4, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
if ( i==3)axis(side=1, las=1)
legend("topright", panlab[i+3], bty="n", cex=1.2) 
#inset=c(-.2,0))


}

#mtext("Turbidity", side = 2, line = 0, outer = TRUE, at = NA,
 #     adj = NA, padj = NA, cex = 1.1, col = NA, font = NA)
mtext("year", side = 1, line = +2, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.4, col = NA, font = NA)

