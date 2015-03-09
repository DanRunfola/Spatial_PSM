#Script to generate Monte Carlo realized dummy datasets and record
#(a) the beta values used to generate the data,  
#(b) the data itself,
#(c) the Moran's I of the data.


library(RandomFields)

#Set a seed so we get the same realization each time..
RFoptions(seed=1602)

#Using an Gaussian model with exponential spatial covariance
model <- RMexp()

#How big the field will be
x <- seq(0, 10, 1)

#Simulate the field (n is number of fields we want.  Here we create 2.)
#These represent an example "outcome", and "control".
z <- RFsimulate(model, x, x, n=2)
dev.off()
plot(z)
z.SPDF <- as(z, 'SpatialPointsDataFrame')

#We also want to simulate a binary treatment, which we do below. 
#Same dimensions, and predicated on the same RMexp model.
#Note this simply takes any values >0 and recodes them as a 1,
#where values are a product of the RMexp() model.
model <- RPbernoulli((model))

#treatment binary.
t <- RFsimulate(model, x, x, n=1)
dev.off()
plot(t)
t.SPDF <- as(t, 'SpatialPointsDataFrame')

f.SPDF <- merge(t,z)
