# In these examples, the data set is 'HD2' and response is 'Anisakas'

########################################################################
# Load packages from R and support functions that we wrote
library(lattice)  
library(mgcv) 
library(MASS)
library(plyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
source("HighstatLibV13.R")
########################################################################


# Zero inflation?
sum(HD2$Anisakis == 0)
sum(HD2$Anisakis == 0) / nrow(HD2)

# The following procedures should be done after running your models to judge how they perform...

# Calculate the dispersion statistic.
E1 <- resid(M1, type = "pearson")
NMinuspa <- M1$df.res 
sum(E1^2) / NMinuspa

#########################################
# Why do we have overdispersion?
# A. Outliers?                  ==> Remove them?
# B. Missing covariates or interactions?  ==> Go back or..add a 
#                                            latent variable 
# C. Zero inflation?            ==> ZIP
# D. Large variance?            ==> NB
# E. Correlation?               ==> GLMM
# F. Non-linear patterns        ==> GAM(M) 
# G. Wrong link function        ==> Change it 


# Model validation M1 
# Plot residuals versus fitted values
par(mfrow = c(1,1))
F1a <- fitted(M1)            
plot(x = F1,                   
     y = E1,                   
     xlab = "Fitted values M1",      
     ylab = "Pearson residuals M1")  
abline(h = 0, lty = 2)            



# Plot residuals vs every covariate (substitute in all your covariates)      
MyVar <- c("TL", "TW", "Month")
HD3$E1 <- E1
MyMultipanel.ggp2(Z = HD3, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

# Can model M1 cope with zero-inflation?
# This is how you can simulate data from a GAM 
simulate(M1)

# Do this 1000 times:
NSim <- 1000
YSim <- matrix(nrow = nrow(HD3), ncol = NSim)
for (i in 1:1000){ 
  YSim[,i] <- simulate(M1a)$sim_1
}

# Calculate the number of zeros in each of the 1,000 data sets.
Zeros <- vector(length = NSim)
for(i in 1:NSim){
  Zeros[i] <- sum(YSim[,i]==0)
}

# Plot this as a table..and do some fancy xlim stuff
RangeTable <- range(as.numeric(names(table(Zeros))))
ZerosInData <- sum(HD3$Anisakis == 0)
par(mfrow = c(1,1))
plot(table(Zeros), 
     xlab = "Number of zeros",
     ylab = "Frequency",
     xlim = c(min(RangeTable[1], ZerosInData), max(RangeTable[2], ZerosInData)))
points(x = sum(HD3$Anisakis == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
# The red dot is the number of zeros in the original data set.
# The frequency plot (black lines) shows the number of zeros in simulated
# data from the Poisson GAM
# Conclusion: the Poisson GAM M1a cannot cope with the 25% of zeros.


# Plot the smoothers 
par(mfrow = c(1,1))
plot(M1)

# We may want to:
#   -reduce the amount of smoothing a little bit
#   -drop any outliers..., etc

# Options to continue:
#  1. ZIP GAM
#  2. NB GAM


# Following codis borrowed from another exercise using another data set 'TT' 
# Check spatial dependency with a sample-variogram
mydata <- data.frame(E2, TT$Xkm, TT$Ykm)
coordinates(mydata) <- c("TT.Xkm", "TT.Ykm")
V12 <- variogram(E2 ~ 1, cutoff = 200, mydata, cressie = TRUE)
plot(V12, 
     main = "Correlation", 
     xlab = "Distance", 
     ylab = "Semi-variogram")

# If you see the points follow a curve from low to high and then level off, that means you have spatial correlation



# Can we visualize that?
MyCol <- c(1,2)[as.numeric(E2>=0) + 1]
MyCex <- 4 * abs(E2) / max(abs(E2)) 
xyplot(Ykm ~ Xkm, 
       aspect = "iso", 
       col = MyCol,
       cex = MyCex,
       xlab = "X coordinate",
       ylab = "Y coordinate",
       data = TT,
       pch = 1)
# Difficult to see...but are the big red dots clustered?
# Based on the variogram we would say that we need a model with spatial 
# correlation, and perhaps also 1 smoothing function for NDVI.



# For this next part, you can view the random effects on a map of the sites to see if there are any patterns. In this code, I3 was the name of the model and 'Nest' was the random effect
# Can we also see this pattern in the random effects of the model I3?
ai    <- I3$summary.random$Nest[,"mean"]
MyCex <- 5 * abs(ai) / max(ai)
MyCol <- ifelse(sign(ai) == 1, 2, 5)
Ykm   <- tapply(Owls$Ykm, FUN = mean, INDEX = Owls$Nest) 
Xkm   <- tapply(Owls$Xkm, FUN = mean, INDEX = Owls$Nest)  
xyplot(Ykm ~ Xkm,
       aspect = "iso",
       cex = MyCex,
       col = MyCol,
       pch = 16)

