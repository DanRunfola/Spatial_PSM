#Script to generate Monte Carlo realized dummy datasets and record
#(a) the beta values used to generate the data,  
#(b) the data itself,
#(c) the Moran's I of the data.

#Small notes for Ariel:
#View(f.SPDF@data)  - allows you to view (capital V) the attributes table within a shapefile, where
#f.SPDF is the shapefile

#install.packages("RandomFields") - how to install a library you're missing.  Only done once, just like
#a module in STATA.

library(RandomFields)
library(sp)
library(spdep)
library(psych)

#User Settings
#Generate a random field with spatial effects in Covariates 
#if 1, covariates have spatial effects.  If 0, completely random fields with no spatial effects.
cov_spatial_effects = 0

#How big the field will be
x <- seq(1, 10, 1)

#--------------------------------
#Initialize Spatial Field
#--------------------------------
#Set a seed if desired (NA for random each iteration)
RFoptions(seed=NA)

#Using an Gaussian model with exponential spatial covariance
model <- RMexp()

#--------------------------------
#Generate Covariates with Spatial Effects
#--------------------------------
if(cov_spatial_effects == 1)
{
  #Simulate the field (n is number of fields we want.  Here we create 2.)
  #These represent an example "outcome", and "control".
  z <- RFsimulate(model, x, x, n=4)
  z.SPDF <- as(z, 'SpatialPolygonsDataFrame')
  
  t <- RFsimulate(model, x, x, n=1)
  t.SPDF <- as(t, 'SpatialPolygonsDataFrame')
  
  f.SPDF = z.SPDF
  f.SPDF@data[5:7]=data.frame(t.SPDF)
  
}

if(cov_spatial_effects == 0)
{
  #Initialize field with an arbitrary run
  z <- RFsimulate(model, x, x, n=1)
  z.SPDF <- as(z, 'SpatialPolygonsDataFrame')
  
  #For each variable, randomly define and name.
  j = 1
  while (j != 6)
  {
    z.SPDF@data[j] <- sample(10000, size = nrow(z.SPDF@data), replace = TRUE)
    j = j + 1
  }
  
  
  f.SPDF <- z.SPDF


}

#Rename the Variables for Later Interpretation
colnames(f.SPDF@data)[1] <- "ControlA"
colnames(f.SPDF@data)[2] <- "RandomFieldA"
colnames(f.SPDF@data)[3] <- "ControlB"
colnames(f.SPDF@data)[4] <- "RandomFieldB"
colnames(f.SPDF@data)[5] <- "Treatment"

#Redefine our Treatment as contingent upon a control and a random field.


f.SPDF@data["Treatment"] = f.SPDF@data["ControlA"] + f.SPDF@data["RandomFieldA"]

#Scale data to calculate treatment binary.
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))

#Convert the Treatment to a Binary
f.SPDF@data$Treatment[which(f.SPDF@data$Treatment > 0)] <- 1 
f.SPDF@data$Treatment[which(f.SPDF@data$Treatment == 0)] <- 1 
f.SPDF@data$Treatment[which(f.SPDF@data$Treatment < 0)] <- 0 
spplot(f.SPDF[1:5])

#--------------------------------------------------------------------------
#Constructing weights matrices and creating data from SLRM
#--------------------------------------------------------------------------

#Neighbor weights sum to 1 for each unit (style=W) to avoid oddities in edge cases, 
#Neighbor is "queens", or any touching unit at one lag.
f.NB = poly2nb(f.SPDF)
f.W = nb2listw(f.NB, style='W')

#Re-scale data (Treatment Binary to SD; note oddity in interpretation**)
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))

#yiA = intercept + (Theta * Treatment) + (Beta * Control (time 2, for now))
yiAfunc <- function(a, b) (0.0+(1.0*a)+(1.0*b))
f.SPDF@data["yiA"] = apply(f.SPDF@data[,c('Treatment','ControlA')], 1, function(y) yiAfunc(y['Treatment'],y['ControlA']))

#Test the Moran's I.  At this stage, it should always be close to 0 and not significant.
moran.test(f.SPDF@data[,6], f.W)

#Re-standardize the new yiA
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))

#Calculate the lagged y variable (from our initial estimation of yiA)
yiA_lag <- lag.listw(f.W, f.SPDF@data$yiA)
f.SPDF@data[7]=data.frame(yiA_lag)

#Now, we add our yiA to the weighted yiA lag term.
#the rho value for the lag term is randomized to test for 
#different degrees of spatial correlation.
#Uniform sampling between 0.1 and 10 for rho.
rho = runif(1,0.1,10)

yiBfunc <- function(a, b) (a + (rho * b))
f.SPDF@data["yiB"] = apply(f.SPDF@data[,c('yiA','yiA_lag')], 1, function(y) yiBfunc(y['yiA'],y['yiA_lag']))

#Re-standardize the new yiB
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))
spplot(f.SPDF[6:8])

#New Moran's test - for higher values of rho, this should be closer to 1 and significant.
moran.test(f.SPDF@data[,8], f.W)

#--------------------------------------------------------------------------
#PSM
#--------------------------------------------------------------------------

#Fit a simple linear model to predict the treatment based on our control variable.
PSM_model <- lm(Treatment ~ ControlA, f.SPDF@data)

plot(f.SPDF@data$Treatment, f.SPDF@data$ControlA)
abline(PSM_model)

#Need to save fitted PSM results...

#Pre PSM Balance
describeBy(f.SPDF@data, group=f.SPDF@data$Treatment)