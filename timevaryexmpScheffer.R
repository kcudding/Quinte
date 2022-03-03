library(deSolve)


#create artificial timeseries
times <- seq(0, 50, by = 0.1)
signal <- as.data.frame(list(times = times, import = rep(0, length(times))))
signal$import <- ifelse((trunc(signal$times) %% 2 == 0), 0, 1)
signal$import=-0.02*signal$times+1.2
signal[8:12,]


         
#Create the interpolating function, using approxfun        
input <- approxfun(signal, rule = 2)
#input(seq(from = 0.98, to = 1.01, by = 0.005))

#Use interpolation function in ODE function
SPCmod <- function(t, x, parms) {
 with(as.list(c(parms, x)), {
 import <- input(t)
  a <- import
  dx=((r*x^p)/(x^p+h^p)+a-b*x)
  res <- c(dx)
  list(res, signal = import)
  })
}

#parameters
a= 0.64254962
#a=.779
p=12
h=1.0
r=.6
b=1.0
sigma=0.05
start=1.7
#start=0.83


parms <- c(a = a, p =p, h = h, r = r, b = b)
xstart <- c(x=start)
out <- ode(y = xstart, times = times, func = SPCmod, parms)

lout=as.data.frame(out)
lout$year=lout$time+1972
plot(y=lout$x, x=lout$year)