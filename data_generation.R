#Script to generate Monte Carlo realized dummy datasets and record
#(a) the beta values used to generate the data,  
#(b) the data itself,
#(c) the Moran's I of the data.


library(RandomFields)
require(sp)

#Set a seed so we get the same realization each time..
RFoptions(seed=1602)

#Using an Gaussian model with exponential spatial covariance
model <- RMexp()

#How big the field will be
x <- seq(1, 10, 1)

#Simulate the field (n is number of fields we want.  Here we create 2.)
#These represent an example "outcome", and "control".
z <- RFsimulate(model, x, x, n=4)
dev.off()
plot(z)
z.SPDF <- as(z, 'SpatialPolygonsDataFrame')

#We also want to simulate a binary treatment, which we do below. 
#Same dimensions, and predicated on the same RMexp model.
#Note this simply takes any values >0 and recodes them as a 1,
#where values are a product of the RMexp() model.
model <- RPbernoulli((model))

#treatment binary.
t <- RFsimulate(model, x, x, n=1)
dev.off()
plot(t)
t.SPDF <- as(t, 'SpatialPolygonsDataFrame')

f.SPDF = z.SPDF
f.SPDF@data[5:7]=data.frame(t.SPDF)

#Rename the Variables for Later Interpretation
colnames(f.SPDF@data)[1] <- "Outcome_T1"
colnames(f.SPDF@data)[2] <- "Outcome_T2"
colnames(f.SPDF@data)[3] <- "Control_T1"
colnames(f.SPDF@data)[4] <- "Control_T2"
colnames(f.SPDF@data)[5] <- "Treatment_Binary"

spplot(f.SPDF[1:5])


