#Script to generate Monte Carlo realized dummy datasets and record
#(a) the beta values used to generate the data,  
#(b) the data itself,
#(c) the Moran's I of the data.


library(sp)
library(spdep)

#Generate an arbitrary 10x10 grid.
gt <- GridTopology(c(0.5, 0.5), c(1, 1), c(10,10))
SP <- as(as(SpatialGrid(gt), "SpatialPixels"), "SpatialPolygons")
class(SP)
plot(SP, axes=TRUE) 

#Define the weights matrix - simple queens (touching neighbors)
#Results in 684 nonzero links (out of 10,000 possible)
ex.nb <- poly2nb(SP, queen=TRUE)
#Note the style=W standardizes weights for each unit predicated on the number of
#neighbors, so weights always sum to 1 for an observation.  "B" allows for binary 1/0 instead.
ex.W <- nb2listw(ex.nb, style="W")

#Data Generation Process will eventually be here..
# set.seed(1)
# SP10 <- SP[sample(length(slot(SP, "polygons")), 10)]
# class(SP10)
# plot(SP10, col="green", add=TRUE) 