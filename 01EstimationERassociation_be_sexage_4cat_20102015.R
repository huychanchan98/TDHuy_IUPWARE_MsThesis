################################################################################
# SUPPLEMENTAL MATERIAL of the article:
#   "A hands-on tutorial on a modelling framework for projections of climate 
#    change impacts on health"
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini.
#
# This code reproduces the analysis described in the tutorial.
# R Code and data released under GNU General Public License (version 3).
#
# [date]
# * an updated version of this code compatible with future
#   versions of the software is available at the personal website of the
#   last author (https://github.com/gasparrini/....)
#
# Requirements:
#  - R version 3.5.1
#  - 'dlnm' package version 2.3.5
#  - 'splines' package version  3.5.1
#  - 'MASS' package version 7.3-50
#
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(splines) ; library(MASS)

################################################################################
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
################################################################################

# Use the observed daily temperature-mortality series to estimate
#     the coefficients defining the exposure-response association.
# In this example, we estimate the temperature-related mortality association 
#     using data from London (UK) between 1990 and 2012.
# ADAPTATION USE DATA FROM BELGIUM 2010-2015

# LOAD OBSERVED DATA - DAILY TEMPERATURE-SERIES BETWEEN 1990-2012 IN LONDON
#  NB: file "lndn_obs.csv" provided as downloadable supplemental material  
#      along this code. Description of the variables available in 
#     "VarDescr_lndn_obs.csv".
setwd("D:/Belgium/OneDrive - KU Leuven/KUL/Thesis/code/tutorial/be2020")

obs <- read.csv("bel_obs_age_2010-2015.csv",colClasses=c(date="Date")) 
#obs <- read.csv("bel_obs_female_age_2010-2015.csv",colClasses=c(date="Date")) 
#obs <- read.csv("bel_obs_female_age_2010-2015.csv",colClasses=c(date="Date")) 

# "bel_obs_age_2010-2015.csv" = age group over Belgium 2010-2015
# "bel_obs_male_age_2010-2015.csv" = sexage (MALE)  group over Belgium 2010-2015
# "bel_obs_female_age_2010-2015.csv" = sexage (FEMALE)  group over Belgium 2010-2015

# DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
# - SPECIFICATION PARAMETERS OF THE EXPOSURE-RESPONSE DIMENSION OF THE CROSS-BASIS

# argvar: main model, cubic natural spline with three internal knots in 
#   the 10th, 75th, 90th percentiles of the temperature distribution

argvar <- list(fun="ns", knots = quantile(obs$tmean,c(10,75,90)/100, na.rm=T),
               Bound=range(obs$tmean,na.rm=T))
# argavar1 & argvar2: additional tests, argavar1 fits a linear function, and
#   argvar2 fits a double threshold, with linear relationship for temperatures 
#   below and above the 10th and 90th percentiles,
#   with null association in between.

argvar1 <- list(fun="lin")
argvar2 <- list(fun="thr", thr.value=quantile(obs$tmean,c(10,90)/100, na.rm=T),
                side="d")
# - SPECIFICATION PARAMETERS OF THE LAG-ASSOCIATION DIMENSION OF THE CROSS-BASIS
# Definition of the maximum lag, that is, 7 days
maxlag <- 7
# arglag: main model, it fits a cubic natural spline with three internal knots 
#   equally-spaced in the log-scale.
arglag <- list(fun="ns",knots=logknots(maxlag,nk=3))

# - CREATE CROSSBASIS OBJECTS
cb <- crossbasis(obs$tmean,maxlag,argvar,arglag)  
cb1 <- crossbasis(obs$tmean,maxlag,argvar1,arglag)   
cb2 <- crossbasis(obs$tmean,maxlag,argvar2,arglag)
# FIT THE MODEL
# Include in the model the crossbasis term of temperature, along with the 
#    indicator for day of day of the week (dow) and natural cubic spline of 
#    time with 8 df per year.

m <- glm(all ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
         data=obs, family=quasipoisson)
m1 <- glm(all ~ cb1 + dow + ns(date,df=round(8*length(date)/365.25)), 
          data=obs, family=quasipoisson)
m2 <- glm(all ~ cb2 + dow + ns(date,df=round(8*length(date)/365.25)), 
          data=obs, family=quasipoisson)
# FIT THE MODEL for different age categories
# Apply the same model defined in the script 01 to the age-specific counts:

m0_64 <- glm(all_0_64 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
             data=obs, family=quasipoisson)
m65_74 <- glm(all_65_74 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
              data=obs, family=quasipoisson)
m75_84 <- glm(all_75_84 ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
              data=obs, family=quasipoisson)
m85plus <- glm(all_85plus ~ cb + dow + ns(date,df=round(8*length(date)/365.25)), 
               data=obs, family=quasipoisson)

# GET PREDICTIONS & PLOT
# - DEFINE PROVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
#    NB: 'varcen' is initially set to a provisional temperature value within 
#        the observed range. 
varcen <- 21 

# - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION
#     FOR EACH AGE CATEGORY

cp0_64 <- crosspred(cb,m0_64,cen=varcen,by=0.1)
cen0_64 <- cp0_64$predvar[which.min(cp0_64$allRRfit)] 

cp65_74 <- crosspred(cb,m65_74,cen=varcen,by=0.1)
cen65_74 <- cp65_74$predvar[which.min(cp65_74$allRRfit)] 

cp75_84 <- crosspred(cb,m75_84,cen=varcen,by=0.1)
cen75_84 <- cp75_84$predvar[which.min(cp75_84$allRRfit)] 

cp85plus <- crosspred(cb,m85plus,cen=varcen,by=0.1)
cen85plus <- cp85plus$predvar[which.min(cp85plus$allRRfit)] 

# - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
pred0_64 <- crosspred(cb, m0_64, cen=cen0_64, by=1)   
pred65_74 <- crosspred(cb, m65_74, cen=cen65_74, by=1)   
pred75_84 <- crosspred(cb, m75_84, cen=cen75_84, by=1)  
pred85plus <- crosspred(cb, m85plus, cen=cen85plus, by=1)  


# - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION 
# MMT corresponds to the temperature of minimum mortality, which will be used as
#    as reference to estimate relative risks and as temperature threshold 
#    to differentiate the contribution of heat and cold to the total mortality 
#    attributable to non-optimal temperatures.

cp <- crosspred(cb,m,cen=varcen,by=0.1)
cen <- cp$predvar[which.min(cp$allRRfit)] 

# - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
pred <- crosspred(cb, m, cen=cen, by=1)   

pred1 <- crosspred(cb1, m1, cen=cen, by=0.1)   
pred2 <- crosspred(cb2, m2, cen=cen, by=0.1)   

################################################################################
# PLOT - FIGURE 1
#pdf("Fig1_bel_age_2010-2015.pdf",height=4.5,width=10)
#pdf("Fig1_bel_male_age_2010-2015.pdf",height=4.5,width=10)
#pdf("Fig1_bel_female_age_2010-2015.pdf",height=4.5,width=10)

xlab <- expression(paste("Temperature (",degree,"C)"))

layout(matrix(c(1,2,3),ncol=3,byrow=T))

# PLOT - 3D
par(mar=c(2,3,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"3d",ltheta=150,xlab="Temperature (C)",ylab="Lag",zlab="RR", 
     col=gray(0.9), main="Exposure-lag-response")

# OVERALL
# The plots show the cumulative exposure-response association, in terms of 
#    relative risks (RR) and centered in the MMT, across the 7 days of lag.
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",main="Overall")

# OVERALL DIFFERENT FUNCTIONS
# See the different shapes of the exposure-response association using the
#   three functions (non-linear, linear, double threshold).
par(mar=c(5,4,3,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",ci="n",main="Overall - different functions")
lines(pred1, col="blue")
lines(pred2, col="forestgreen")

legend("top",c("Natural spline (Main)","Linear","Double threshold")
       ,xpd = TRUE,col=c("red","blue","forestgreen"),lwd=2,bg="white"
       ,cex=0.8,ncol=1,inset=0.04)

layout(1)
#dev.off()

# NB: run pdf() and dev.off() for saving the plot in pdf format.
#############################################################################
#FIGURE 6
#pdf("Fig6_bel_age_2010-2015.pdf",height=5,width=9)
#pdf("Fig6_bel_male_age_2010-2015.pdf",height=5,width=9)
#pdf("Fig6_bel_female_age_2010-2015.pdf",height=5,width=9)

# FIGURE 6A DEMOGRAPHIC CHANGES
layout(matrix(c(1,2),ncol=2,byrow=T))
xlab <- expression(paste("Temperature (",degree,"C)"))

par(mar=c(5,4,4,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)
plot(pred0_64,"overall",col="magenta",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",ci="n",main="Demographic changes", lwd=2)
lines(pred65_74, col="blue", lwd=2, lty=2)
lines(pred75_84, col="forestgreen", lwd=2, lty=3)
lines(pred85plus, col="orange", lwd=2, lty=4)

legend(17,2.4,c("<65 years","65-74 years","75-84 years", ">84 years")
       ,xpd = TRUE,col=c("magenta","blue","forestgreen", "orange"),lwd=2,bg="white"
       ,cex=0.9,ncol=1,inset=0.04,bty="n", lty=1:4)
text(21,2.5, "Age-specific relationships", cex=1)


# FIGURE 6B ADAPTATION

# INDICATOR FOR HEAT (ABOVE MMT)
ind <- pred$predvar>=cen 
par(mar=c(5,4,4,1),mgp=c(3,1,0),las=1,cex.axis=0.9,cex.lab=1)

plot(pred,"overall",col="red",ylim=c(0.5,2.5),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",ci="n",main="Adaptation")

lines(pred$predvar,pred$allRRfit,col=2,lwd=2)
# PLOT LINE OF THE MODIFIED CURVE (30% REDUCTION OF HEAT) 
lines(pred$predvar[ind],exp(pred$allfit[ind]/1.7),col=2,lwd=2,lty=2)
legend(17,2.4,c("No adaptation","With adaptation"),lty=1:2,inset=0.005,
       cex=0.9,bty="n",col=2)
text(21,2.5, "Reduction in heat risk", cex=1)

layout(1)
#dev.off()


################################################################################





