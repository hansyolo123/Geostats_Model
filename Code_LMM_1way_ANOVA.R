##############################################
##### Fitting LMM w/ correlated error ########
##############################################

#Focus is 1 way ANOVA w/correlated errors

setwd("/Users/tfabishop/Dropbox/Research/Code/LMM_Univariate/R/Geostats_Model")

#Read data
data<-read.csv("~/Dropbox/Research/Code/LMM_Univariate/Data/data.csv")

head(data)
str(data)

summary(data$DSM1)

tapply(data$ph_0_30,data$DSM1,summary)

install.packages("gstat")

library(gstat)

library(sp)

#Set up model
base.lm<-lm(ph_0_30~DSM1,data=data)

#DIAGNOSTICS
full<-base.lm   #tell R the object where the regression model is stored
par(mfrow=c(2,2))
hist(rstandard(full),xlab="Standardised residuals",ylab="Frequency",main=NULL)
plot(full$fitted.values,(rstandard(full)),xlab="Fitted values",
     ylab="Standardised residuals")
qqnorm(rstandard(full),main=NULL,xlab="Normal quantiles",ylab="Sample quantiles")
abline(0,1)

names(clay.lm)

data$residuals<-base.lm$residuals

soils<-SpatialPointsDataFrame(data[,4:5],data)

summary(soils)

# Object of class SpatialPointsDataFrame
# Coordinates:
#   min     max
# east   376490  403486
# north 6658420 6684736
# Is projected: NA 
# proj4string : [NA]
# Number of points: 129
# Data attributes:
#   east            north              clay            em34          residuals      
# Min.   :376490   Min.   :6658420   Min.   : 5.30   Min.   : 55.23   Min.   :-45.502  
# 1st Qu.:383801   1st Qu.:6667801   1st Qu.:41.06   1st Qu.:101.49   1st Qu.: -5.861  
# Median :388538   Median :6672662   Median :50.94   Median :140.16   Median :  2.080  
# Mean   :388670   Mean   :6671760   Mean   :46.86   Mean   :137.93   Mean   :  0.000  
# 3rd Qu.:392778   3rd Qu.:6675810   3rd Qu.:57.98   3rd Qu.:170.49   3rd Qu.:  7.138  
# Max.   :403486   Max.   :6684736   Max.   :70.15   Max.   :211.66   Max.   : 36.909  

bubble(soils,"residuals")

###Variogram Modelling

resid.vcld <- variogram(residuals ~ 1, soils, cloud=TRUE)

plot(variogram(residuals ~ 1, soils, cloud=TRUE))

resid.vgm <- variogram(residuals ~ 1, soils)

plot(resid.vgm)

resid.sph <- fit.variogram(resid.vgm,model=vgm(1,"Sph",10000,1))

resid.exp <- fit.variogram(resid.vgm,model=vgm(1,"Exp",10000,1))

attr(resid.sph, 'SSErr')

#[1] 0.06596707

attr(resid.exp, 'SSErr')

#[1] 0.07249274

#Spherical has lower RSS so we choose that model

plot(resid.vgm,resid.sph)

resid.sph

#   model    psill    range
# 1   Nug 83.88134    0.000
# 2   Sph 84.27120 9046.839


#####Part 3 - Linear mixed models in R

  

install.packages("geoR")

install.packages("raster")

install.packages("epiR")

library(geoR)
library(lattice)
library(gstat)
library(moments)
library(ggplot2)
library(raster)
library(epiR)

citation()###Exploratory Data Analysis

soil.geodata <- as.geodata(soil,coords.col=4:5,data.col=30, covar.col=22)

summary(soil.geodata)

# $n
# [1] 129
# 
# $coords.summary
# x        y
# min 0.027411 0.100981
# max 2.727011 2.732581
# 
# $distances.summary
# min         max 
# 0.006946222 3.005036456 
# 
# $data.summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.30   41.06   50.94   46.86   57.98   70.15 
# 
# $covariate.summary
# clay            em34              x                 y        
# Min.   : 5.30   Min.   : 55.23   Min.   :0.02741   Min.   :0.101  
# 1st Qu.:41.06   1st Qu.:101.49   1st Qu.:0.75851   1st Qu.:1.039  
# Median :50.94   Median :140.16   Median :1.23221   Median :1.525  
# Mean   :46.86   Mean   :137.93   Mean   :1.24538   Mean   :1.435  
# 3rd Qu.:57.98   3rd Qu.:170.49   3rd Qu.:1.65621   3rd Qu.:1.840  
# Max.   :70.15   Max.   :211.66   Max.   :2.72701   Max.   :2.733  

plot(soil.geodata)

#Experimental semivariogram

soil.var <- variog(soil.geodata, trend = ~DSM1)

par(mfrow=c(1,1))

plot(soil.var,pch = 19, col= "blue", main="",xlim=c(0,1),ylim=c(0,30000),xlab="distance(m)")

#Initial parameter estimates

ini.=matrix(NA,ncol=2,nrow=9)
ini.[,1]=c(1000,5000,20000,1000,5000,20000,1000,5000,20000)
ini.[,2]=c(0.5,0.5,0.5,1.0,1.0,1.0,1.5,1.5,1.5)

lmm.exp <- likfit(soil.geodata, trend= ~DSM1,lambda=1,ini.cov.pars = ini.,fix.nugget = FALSE, nugget = 10, 
                  lik.method = "REML", cov.model = "exp")

lmm.sph <- likfit(soil.geodata, trend= ~DSM1,lambda=1,ini.cov.pars = ini.,fix.nugget = FALSE, nugget = 10, 
                  lik.method = "REML", cov.model = "sph")

summary(lmm.exp)

# Maximised Likelihood:
#   log.L n.params      AIC      BIC 
# "-484.9"      "7"  "983.8"   "1004" 
# 
# non spatial model:
#   log.L n.params      AIC      BIC 
# "-492.3"      "5"  "994.6"   "1009" 

summary(lmm.sph)

# Maximised Likelihood:
#   log.L n.params      AIC      BIC 
# "-484.4"      "7"  "982.8"   "1003" 
# 
# non spatial model:
#   log.L n.params      AIC      BIC 
# "-492.3"      "5"  "994.6"   "1009"

#fitted variogram

plot(soil.var,pch = 19, col= "blue", main="",xlim=c(0,1),ylim=c(0,25000),xlab="distance(m)")
lines(lmm.exp, col="red")
lines(lmm.sph, col="green")

#For 2 df difference the chisq critical value
qchisq(0.95, df=2)

#Spherical model best

#Spatial vs Non-Spatial
pchisq(q=-484.4--492.3,df=2,lower.tail=FALSE)

# [1] 0.0192547

###Variable selection

fitmodel<-lmm.sph
datas<-data

##Code for getting P-values from geoR object
## get coefficients
coefficients <- fitmodel$beta
## get standard errors
se_error <- sqrt(diag(fitmodel$beta.var))
## get t values
t_value <- coefficients/se_error
## and probabilities
t_prob <- 2 * pt(-abs(t_value), df = (nrow(datas) -length (coefficients)))
## make pretty
coef_mat <- cbind(coefficients, se_error, t_value, t_prob)
colnames(coef_mat) <- c("Estimate", "Std.Err","t value", "Pr(>|t|)")
rownames(coef_mat) <- colnames(model.matrix(fitmodel$trend, soil.geodata))
printCoefmat(coef_mat)

#                 Estimate   Std.Err t value  Pr(>|t|)    
#   (Intercept) 20.982494  7.731520  2.7139  0.007589 ** 
#   x            1.974214  2.676611  0.7376  0.462152    
#   y           -4.943037  2.754157 -1.7948  0.075109 .  
#   em34         0.213832  0.035942  5.9493 2.524e-08 ***

#Compare to OLS model

summary(base.lm)

#   Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 18.48838    5.52570   3.346  0.00108 ** 
#   x            2.08798    1.81041   1.153  0.25098    
#   y           -3.97289    1.92893  -2.060  0.04151 *  
#   em34         0.22818    0.02784   8.197  2.5e-13 ***

