#Exploring the Moran Eigenvector 

#Very simple case - 3 units of observation, A B and C
#A borders B, B borders A and C, C borders A.
# [A][B][C]

#First define the Weights Matrix in the simplest way (contiguity)
cont_vector <- c(0,1,2,1,0,1,2,1,0)
rownames_vector = c("A","B","C")
colnames_vector = c("A","B","C")
W <- matrix(cont_vector,nrow=3,ncol=3,byrow=TRUE, dimnames=list(rownames_vector,colnames_vector ))
W

#Define an Identity Matrix (could also use diag(3))
cont_vector <- c(1,0,0,0,1,0,0,0,1)
rownames_vector = c("A","B","C")
colnames_vector = c("A","B","C")
I <- matrix(cont_vector,nrow=3,ncol=3,byrow=TRUE, dimnames=list(rownames_vector,colnames_vector ))
I

#Define an X vector 
X <- matrix(c(1,1,1), nrow=3, ncol=1)
X

#These are equivalent:
#Specification in Roger Bivand, R users guide
L = I - X%*%solve(t(X)%*%X)%*%t(X)
M = L %*% W %*% L
M

#Specification in Julie Le Gallo and Antonio Paex (Environment and Planning A)
M = (I - X%*%t(X)/3) %*%W%*% (I - X%*%t(X)/3)
M

#NOTE ON THE ABOVE: It is not clear if the matrices should be multiplied using element-wise multiplication
#or, matrix multiplication.  It is assumed the intent was matrix multiplication, which is likely a good guess.
#However, the interpretation of the eigenvectors is not as detailed in Gallo/Paez in that case,
#which indicates something may be amiss.

#The next step in this is thus to grab the R program and compare the outputs of this very simple script,
#and see how they do it (SpatialFiltering)

#Ultimately, through either calculation M is a projection matrix.
eigen(M)
