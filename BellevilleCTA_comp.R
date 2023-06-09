library(strucchange)
library(ecotraj)

rm(list=ls())
#read in data supplied by DFO, rename and aggregate into annuual means
library(openxlsx)
library("readxl")
my_data <- read.xlsx("F:/homeoffice/kim/Quinte/QUINTE.xlsx", sheet = 1)
my_data2 <- read.xlsx("F:/homeoffice/kim/Quinte/All phyto annual average Comp and Ind date avgwith count.xlsx", 
                      sheet=3, startRow=4)
#rename col in recent file
colnames(my_data)[colnames(my_data)=="Epar#Ship"] <- "epar"
colnames(my_data)[colnames(my_data)=="Secchi_depth#Ship"] <- "secchi"

phytos=c("Cyanophyta",  "Chlorophyta",
         "Euglenophyto", "Chrysophyceae","Diatomeae" ,
         "Cryptophyta" , "Dinophyceae" )
#[65] "Cyanophyta.biomass"                             
#[66] "Chlorophyta.biomass"                            
#[67] "Euglenophyto.biomass"                           
#[68] "Chrysophyceae.biomass"                          
#[69] "Diatomeae.biomass"                              
#[70] "Cryptophyta.biomass"                            
#[71] "Dinophyceae.biomass" 


#mmetric="Chlorophyta.biomass" #select which metric to analyze

lobs<-function(x) sum(!is.na(x))

i <- c(65:71)  
my_data[ ,i] <- apply(my_data[ ,i], 2,function(x) as.numeric(as.character(x)))
colnames(my_data)[colnames(my_data)=="Chl_a_uncorrected#GLLFAS"]<-"Chla"

  agg <- aggregate(list(metric=my_data[, c(15,52,65:71)],TP=my_data$TP), 
                   by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                   FUN=mean, na.rm=TRUE, na.action="na.pass")
  colnames(agg)[3:11]=paste0(c("Chla", "epar",phytos),".mean")
  agg2 <- aggregate(list(metricN=my_data[,c(15,52,65:71)]), 
                    by = list(Station_Acronym=my_data$Station_Acronym, year=my_data$year), 
                    FUN=lobs)
  colnames(agg2)[3:11]=paste0(c("Chla", "epar",phytos),".n")
  
  trym=merge(agg,agg2)
  
  aggwgts <- apply(trym[,13:21], 2, function(x) x/mean(x))
  colnames(aggwgts)=paste0(c("Chla", "epar",phytos),".wgt")
  comm=cbind(trym, aggwgts)
  
  
  #select sites on which to run analysis
  stvec=c("B")
  
  stlabel=c("Belleville")
  #select sites on which to run analysis
  
  
  spydat=comm[comm$Station_Acronym%in%stvec,]
  cdat=my_data2[my_data2$station%in%stvec,]
  
  xyear=c(1980,2011)
  spydat=spydat[spydat$year%in%c(1980:2011),]
  spydat$Station_Acronym=as.factor(spydat$Station_Acronym)
  spyyear=spydat[complete.cases(spydat),]
  colnames(spyyear)[1]="station"
 # spyyear=cdat
  rownames(spyyear)=spyyear$year
  randspy=spyyear
  randspy$station="R"
  randspy[,3:9]=apply(spyyear[,3:9], 2, function(x) sample(x))
 comm=rbind(spyyear, randspy)
  D = dist(spyyear[-6,c(3:4, 6:9)])
  rm87=which(spyyear$year==1987)
  D2 = dist(log(spyyear[-rm87,3:9]+1))
  D = dist(log(spyyear[-6,c(3:4, 6:9)]+1))
 # D2 = dist((spyyear[,3:9]+1))
  cmdscale(D)
  cmd_D2=cmdscale(D, eig = TRUE, add = TRUE, k = nrow(as.matrix(D))-1)
  x <- cmd_D2$points[,1]
  y <- cmd_D2$points[,2]
  spyyear$col=ifelse(spyyear$year>1994, 2, 1)
  plot(x, y, type = "p", pch=16, col=spyyear$col,
       xlab = paste0("PCoA ", 1,
            " (", round(100 * cmd_D2$eig[1]/sum(cmd_D2$eig)), "%)"),
  ylab = paste0("PCoA ", 2, 
      " (", round(100 * cmd_D2$eig[2]/sum(cmd_D2$eig)), "%)"))

  text(x,y, labels=spyyear$year[-6], cex.lab=0.4, col=spyyear$col[-6],pos=3)
 # , asp = 1,
  xlab = paste0("PCoA ", 1, " (", round(100 * cmd_D2$eig[1]/sum(cmd_D2$eig)), "%)")
#  ylab = paste0("PCoA ", 2, " (", round(100 * cmd_D2$eig[2]/sum(cmd_D2$eig)), "%)"))
  trajectoryPCoA(D, sites=comm$station[-6], surveys=comm$year[-6], 
                 traj.colors=c("black", "red"), lwd = 2,
                 survey.labels = T)
  library(vegan)
  rm87=6
  svarespec = wisconsin(spyyear[-rm87,3:9])
  svarespec=decostand(spyyear[-rm87,3:9], method = "hellinger")
  disimvar = vegdist(svarespec, method = "bray")
  
  PCoA <- cmdscale(disimvar, k = 2, eig = T, add = T)
  x <- PCoA$points[,1]
  y <- PCoA$points[,2]
  
  plot(x, y, type = "p", pch=16, col=spyyear$col,
       xlab = paste0("PCoA ", 1,
                     " (", round(100 * cmd_D2$eig[1]/sum(cmd_D2$eig)), "%)"),
       ylab = paste0("PCoA ", 2, 
                     " (", round(100 * cmd_D2$eig[2]/sum(cmd_D2$eig)), "%)"))
  
  text(x,y, labels=spyyear$year[-rm87], cex.lab=0.4, col=spyyear$col,pos=3)
  rownames(cdat)=cdat$year
  svarespec=decostand(cdat[-rm87,c(3:4,6:9)], method = "hellinger")
  
  ord3 = vegan::rda(svarespec ~ spyyear$TP[-rm87]+spyyear$year[-rm87])
  ord3 = vegan::rda(svarespec ~ spydat$TP[-rm87]+spydat$year[-rm87])
  ord3 = vegan::rda(svarespec ~ spydat$TP[-rm87])
  anova(ord3, by = "term", permutations = 199)
  plot(ord3, scaling = "species", main = "Sites scaling (scaling=1)",
       correlation = TRUE)
  plot(ord3, scaling = "sites", main = "Sites scaling (scaling=1)",
       correlation = TRUE)
  
  
  #TP model
  tp1=(lm(x~spyyear$TP[-rm87]+spyyear$year[-rm87], weights=spyyear$weights[-rm87], data=spyyear[-rm87]))
  tp2=(lm(metric~year, weights=weights, data=spydat[-rm87]))
  tp4=(lm(metric~TP+year, weights=weights, data=spydat[-rm87]))
  
  temp_comm=temp_comm[-2,]
  
  svarespec=decostand(temp_comm[,3:9], method = "hellinger")
  ord3 = vegan::rda(svarespec ~ temp_comm$TP+temp_comm$mean)
  anova(ord3, by = "term", permutations = 199)
  plot(ord3, scaling = "sites", main = "Sites scaling (scaling=1)",
       correlation = TRUE)
  
  disimvar = vegdist(temp_comm[,3:9], method = "euclid")
  PCoA <- cmdscale(disimvar, k = 2, eig = T, add = T)
  
  text(x,y, labels=temp_comm$year, cex.lab=0.4, col=temp_comm$col,pos=3)