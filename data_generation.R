#Script to generate Monte Carlo realized dummy datasets and record
#(a) the beta values used to generate the data,  
#(b) the data itself,
#(c) the Moran's I of the data.

#Small notes for Ariel:
#View(f.SPDF@data)  - allows you to view (capital V) the attributes table within a shapefile, where
#f.SPDF is the shapefile

#install.packages("RandomFields") - how to install a library you're missing.  Only done once, just like
#a module in STATA.

library(RandomFields) #Generating Random Spatial Fields
library(sp)    #Spatial data handling
library(spdep) #Spatial Lag Model Fitting
library(psych) #Helps with comparing by groups
library(FNN) #Nearest Neighbor Classifications
library(rgl) #3D Scatterplots.  Note this requires X11 on Ubuntu, see below command:
#sudo apt-get install libX11-dev freeglut3 freeglut3-dev

#User Settings
#Generate a random field with spatial effects in Covariates 
#if 1, covariates have spatial effects.  If 0, completely random fields with no spatial effects.
cov_spatial_effects = 0

#How big the field will be
x <- seq(1, 10, 1)

#How many iterations to perform
total_iterations = 1000

#Do we create maps? 1 = Yes.
maps = 0

#Do we output plots and iteration counts?
verbose = 0

#----------------------------------
#----------------------------------
#Initialize the Data Frame
beta_df <- data.frame(matrix(ncol=3, nrow=total_iterations))
colnames(beta_df)[1] <- "BdifBhat"
colnames(beta_df)[2] <- "Rho"
colnames(beta_df)[3] <- "Morans_I"

beta_df_NL <- data.frame(matrix(ncol=2, nrow=total_iterations))
colnames(beta_df_NL)[1] <- "BdifBhat"
colnames(beta_df_NL)[2] <- "Morans_I"

it_cnt = 1
while (it_cnt < (total_iterations+1))
{
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
  if(maps == 1)
  {
    spplot(f.SPDF[1:5])
  }
  #--------------------------------------------------------------------------
  #Constructing weights matrices and creating data from SLRM
  #--------------------------------------------------------------------------
  
  #Neighbor weights sum to 1 for each unit (style=W) to avoid oddities in edge cases, 
  #Neighbor is "queens", or any touching unit at one lag.
  f.NB = poly2nb(f.SPDF)
  f.W = nb2listw(f.NB, style='W')
  
  #Re-scale data (Treatment Binary to SD; note oddity in interpretation**)
  f.SPDF@data <- data.frame(lapply(f.SPDF@data, function(x) scale(x)))
  
  #yiA = intercept + (Theta * Treatment) + (Beta * ControlA)
  yiAfunc <- function(a, b) (0.0+(1.0*a)+(1.0*b))
  f.SPDF@data["yiA"] = apply(f.SPDF@data[,c('Treatment','ControlA')], 1, function(y) yiAfunc(y['Treatment'],y['ControlA']))
  
  #Test the Moran's I.  At this stage, it should always be close to 0 and not significant.
  mi_NL = moran.test(f.SPDF@data[,6], f.W)
  
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
  
  if(maps == 1)
  {
    spplot(f.SPDF[6:8])
  }
  #New Moran's test - for higher values of rho, this should be closer to 1 and significant.
  mI = moran.test(f.SPDF@data[,8], f.W)
  
  #--------------------------------------------------------------------------
  #PSM
  #--------------------------------------------------------------------------
  
  #Fit a simple linear model to predict the treatment based on our control variable.
  #Note 1: Reconstruct our Treatment as a Binary, for interpretation.
  #Note 2: Standard Error tests are invalid for the PSM.  I don't believe this matters, 
  #As we are interested in the predictive power of ancillary for inclusion in the treatment.
  #Check on this assumption.  Could easily swap for probit or logit.
  
  
  #Convert the Treatment to a Binary 0/1
  f.SPDF@data$Treatment[which(f.SPDF@data$Treatment > 0)] <- 1 
  f.SPDF@data$Treatment[which(f.SPDF@data$Treatment == 0)] <- 1 
  f.SPDF@data$Treatment[which(f.SPDF@data$Treatment < 0)] <- 0 
  
  PSM_model <- lm(Treatment ~ ControlA, f.SPDF@data)
  
  #plot(f.SPDF@data$Treatment, f.SPDF@data$ControlA)
  #abline(PSM_model)
  
  #Need to save fitted PSM results...
  f.SPDF@data["PSM"] <- predict(PSM_model, f.SPDF@data)
  
  #Pre PSM Balance
  describeBy(f.SPDF@data, group=f.SPDF@data$Treatment)
  
  #Run a KNN - for later editing this is a custom-coded of a KNN without replacement.
  #Threshold Value for KNN can be set at any arbitrary distance value,
  #where distance is the PSM fit difference (un-standardized)
  
  #Create a new version of the f.SPDF data frame which records only matched pairs.
  #Further, add a column to keep track of the matches.
  m.SPDF = f.SPDF
  m.SPDF$match = 0
  
  #add an ID column to track matches in the m.SPDF frame, and to remove from f.SPDF  
  m.SPDF$PSM_ID <- seq_len(nrow(m.SPDF))
  f.SPDF$PSM_ID <- seq_len(nrow(f.SPDF))
  
  #add a distance column to track the PSM distance for later analysis
  m.SPDF$PSM_distance <- -1
  
  #add a match pair so we can view what was matched with what.
  m.SPDF$PSM_match_ID <- -1
  
  Treatment_n = length(f.SPDF@data$Treatment[which(f.SPDF@data$Treatment == 1)])
  Control_n = length(f.SPDF@data$Treatment[which(f.SPDF@data$Treatment == 0)])
  cnt = min(Treatment_n, Control_n)
  
  #Loop through all treatment cases to find a match.
  for (j in 1:cnt)
  {
    
    #Run the KNN for all neighbors.  We want to optimize the total distance between PSM (treatment, control) 
    #to be as low as possible.  Thus, we run the full set, choose the lowest distance, drop both pairs, repeat.
    k <- get.knnx(f.SPDF@data$PSM[which(f.SPDF@data$Treatment == 1)], f.SPDF@data$PSM[which(f.SPDF@data$Treatment == 0)], 1)
    #print(k)
    
    #Add the matched treatment and control values to the m.SPDF data frame
    best_m = as.matrix(apply(k$nn.dist, 2, which.min))[1]
   
    #Control PSM ID
    Control_ID = f.SPDF@data$PSM_ID[which(f.SPDF@data$Treatment == 0)][as.matrix(apply(k$nn.dist, 2, which.min))[1]]
  
    
    #Treatment PSM ID
    k_match_id = k[[1]][best_m]
    Treatment_ID = f.SPDF@data$PSM_ID[which(f.SPDF@data$Treatment == 1)][k_match_id]
  
    
    #Add the Treatment ID to the Control Row and Record the distance / matchID between pairs
    m.SPDF@data$match[which(m.SPDF@data$PSM_ID == Control_ID)] = Treatment_ID
    m.SPDF@data$PSM_distance[which(m.SPDF@data$PSM_ID == Control_ID)] = k[[2]][best_m]
    m.SPDF@data$PSM_match_ID[which(m.SPDF@data$PSM_ID == Control_ID)] = j
    
    #Add the Control ID to the Treatment Row and Record the distance / matchID between pairs
    m.SPDF@data$match[which(m.SPDF@data$PSM_ID == Treatment_ID)] = Control_ID
    m.SPDF@data$PSM_distance[which(m.SPDF@data$PSM_ID == Treatment_ID)] = k[[2]][best_m]
    m.SPDF@data$PSM_match_ID[which(m.SPDF@data$PSM_ID == Treatment_ID)] = j
    
    #Drop the paired match out of the f.SPDF matrix 
    f.SPDF@data <- f.SPDF@data[f.SPDF@data$PSM_ID != Treatment_ID ,]
    f.SPDF@data <- f.SPDF@data[f.SPDF@data$PSM_ID != Control_ID ,]
    
  }
  
  #Map the Matches
  sp.label <- function(x, label)
  {
    list("sp.text", coordinates(x), label)
  }
  ISO.sp.label <- function(x)
  {
    sp.label(x, x@data["PSM_match_ID"][[1]])
  }
  make.ISO.sp.label <- function(x)
  {
    do.call("list", ISO.sp.label(x))
  }
  if(maps == 1)
  {
    spplot(m.SPDF['PSM_distance'], sp.layout=make.ISO.sp.label(m.SPDF))
  }
  
  #Calculate the Spatial Euclidean Distance for Each Pair to contrast to PSM NN distance.
  m.SPDF['Euclidean_distance'] <- -1
  for (j in 1:length(m.SPDF))
  {
    #Select the pair
    temp.SPDF <- m.SPDF
    t.PSM_match_ID <- m.SPDF@data["PSM_match_ID"][[1]][j]
    temp.SPDF@data <- m.SPDF@data[m.SPDF@data$PSM_match_ID == t.PSM_match_ID ,]
    
    #Calculate the distance
    d.EUC = dist(temp.SPDF@data, method="euclidean")[[1]]
    m.SPDF@data$Euclidean_distance[[j]] <- d.EUC
    
  }
  
  dist.DF <- m.SPDF@data
  dist.DF <- dist.DF[dist.DF$PSM_distance >0 ,]
  
  title = paste("Moran's I:",round(mI[[3]][[1]],3),"rho:",round(rho,3))
  dist_model <- lm(PSM_distance ~ Euclidean_distance, dist.DF)
  if(verbose == 1)
  {
    plot(dist.DF$Euclidean_distance, dist.DF$PSM_distance, main=title)
    abline(dist_model)
  }
  
  
  
  #------------------------------
  #Balancing based on PSM distance
  #Ludicrously simple for now, need to complicate:
  #Keep matches that are less than mean difference - may want to change this later.
  #------------------------------
  describe(dist.DF$PSM_distance)
  #median.PSM = median(dist.DF$PSM_distance)
  mean.PSM = describe(dist.DF$PSM_distance)[3][[1]]
  
  dist.DF <- dist.DF[dist.DF$PSM_distance < mean.PSM ,]
  
  #Check our new PSM Balance:
  describeBy(dist.DF, group=dist.DF$Treatment)
  
  #Plot it:
  dist_model <- lm(PSM_distance ~ Euclidean_distance, dist.DF)
  if(verbose == 1)
  {
    plot(dist.DF$Euclidean_distance, dist.DF$PSM_distance, main=title)
    abline(dist_model)
  }
  
  
  
  #In the DGN, both of beta's are equal to 1, with an intercept of 0.
  #use a linear model to estimate our beta in the NO lag case:
  post_psm_model_noLag = lm(yiA ~ Treatment + ControlA, dist.DF)
  noLag_Treatment_beta = summary(post_psm_model_noLag)[4][[1]][2]
  noLag_Control_beta = summary(post_psm_model_noLag)[4][[1]][3]
  
  post_psm_model_Lag = lm(yiB ~ Treatment + ControlA, dist.DF)
  Lag_Treatment_beta = summary(post_psm_model_Lag)[4][[1]][2]
  Lag_Control_beta = summary(post_psm_model_Lag)[4][[1]][3]
  
  #final elements to record:
  #Different in Treatment Effect Beta's from lagged DGN
  BdifBhat = 1 - Lag_Treatment_beta
  
  BdifBhat_NL = 1 - noLag_Treatment_beta
  
  #Save values from iteration
  beta_df["BdifBhat"][it_cnt,] = BdifBhat
  beta_df["Rho"][it_cnt,] = rho
  beta_df["Morans_I"][it_cnt,] = mI[[3]][[1]]
  
  #Also save the "non lagged" data, to benchmark.  
  beta_df_NL["BdifBhat"][it_cnt,] = BdifBhat_NL
  beta_df_NL["Morans_I"][it_cnt,] = mi_NL[[3]][[1]]
  
  it_cnt = it_cnt + 1
  if(verbose == 1)
  {
    print("")
    print("-------------")
    print("Begin Iteration:")
    print(it_cnt)
    print("-------------")
  }
}

if (verbose == 1)
{
  plot3d(beta_df)
  
  plot(beta_df$Morans_I, beta_df$BdifBhat)
  points(beta_df_NL$Morans_I, beta_df_NL$BdifBhat, col="red")
}