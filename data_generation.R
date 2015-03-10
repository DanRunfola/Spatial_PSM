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
  
  #If you want to define the treatment as a binary, uncomment this line.
  #model <- RPbernoulli((model))

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
colnames(f.SPDF@data)[1] <- "ControlA_T1"
colnames(f.SPDF@data)[2] <- "ControlA_T2"
colnames(f.SPDF@data)[3] <- "ControlB_T1"
colnames(f.SPDF@data)[4] <- "ControlB_T2"
colnames(f.SPDF@data)[5] <- "Treatment"
spplot(f.SPDF[1:5])

#--------------------------------------------------------------------------
#Raw data generation (initiatlized surfaces) complete.
#Constructing weights matrices and creating data from a SLRM now.

#Biggest concerns: (a) not sure what modules to use for spatio-temporal regressions;
#(b) I do not incorporate temporal autocorrelation into any of my spatial fields, now.
#Calculate Spatial Weights Matrix
#Very basic assumptions here -
#Neighbor weights sum to 1 for each unit (style=W) to avoid oddities in edge cases, 
#Neighbor is "queens", or any touching unit at one lag.
f.NB = poly2nb(f.SPDF)
f.W = nb2listw(f.NB, style='W')
#Scale all data to standard deviations from mean for easy comparison later
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))

#follow a two-step functional form:
#yiB = yiA + (pWyA) 

#yiA = intercept + (Theta * Treatment) + (Beta * Control (time 2, for now))
yiAfunc <- function(a, b) (0.0+(1.0*a)+(0.5*b))
f.SPDF@data["yiA"] = apply(f.SPDF@data[,c('Treatment','ControlA_T2')], 1, function(y) yiAfunc(y['Treatment'],y['ControlA_T2']))

#Test the Moran's I.  At this stage, it should always be close to 0 and not significant.
#This is because the current model is predicated entirely on random data.
moran.test(f.SPDF@data[,6], f.W)

#Re-standardize the new yiA
f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))

#Now, we want to add in spatial effects.
#First, we need to calculate the lagged y variable (from our initial estimation of yiA)
yiA_lag <- lag.listw(f.W, f.SPDF@data$yiA)
f.SPDF@data[7]=data.frame(yiA_lag)

#Now, we add our yiA to the weighted yiA lag term.
#the rho value for the lag term is set equal to 0.5
yiBfunc <- function(a, b) (a + (0.5 * b))

