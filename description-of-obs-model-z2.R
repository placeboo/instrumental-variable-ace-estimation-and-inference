rm(list=ls())
### In the note below we perform the following calculation
### We derive bounds on the observable distribution under the IV model

require("rcdd")

p <- 2  ## dimension of Z
# X and Y are binary

vmat <- matrix(0,nrow=4*(2^p-1),ncol=4*p)  ## Num rows = num of extreme points
										   ## See this as number of ways to 
										   ## Partition X(z)'s into those getting 1 and
										   ## getting 0, times values given to
										   ## Y(X(z)), which are either 4 or 2
										   ## depending on whether {X(z); z=1,..,p}
										   ## contains 2 values or just one
										   ## See Bonet for general expression

										   ## Num cols = one for each p(x,y|z)
## Ordering of cols: p(x0,y0|z1) p(x0,y1|z1) p(x1,y0|z1) ... p(x1,y1 | zp)

## Ordering of rows: 
## First two rows where X(z) takes 0 for all z, with Y(X(z))=0,1
## Then: two rows where X(z) takes 1 for all z, with Y(X(z))=0,1
## Then: 4*(2^p-2) rows where X(z) takes two different values,
##        ordering of rows is given by Y(X(z)=0)=0, Y(X(z)=1)=0
##        ordering of rows is given by Y(X(z)=0)=0, Y(X(z)=1)=1
##        ordering of rows is given by Y(X(z)=0)=1, Y(X(z)=1)=0
##        ordering of rows is given by Y(X(z)=0)=1, Y(X(z)=1)=1

vmat[1,] <- rep(c(1,0,0,0),p)
vmat[2,] <- rep(c(0,1,0,0),p)
vmat[3,] <- rep(c(0,0,1,0),p)
vmat[4,] <- rep(c(0,0,0,1),p)

## {1}=0 {2}=1
vmat[5,] <-  c(c(1,0,0,0), c(0,0,1,0)) #(0,0)
vmat[6,] <-  c(c(1,0,0,0), c(0,0,0,1)) #(0,1)
vmat[7,] <-  c(c(0,1,0,0), c(0,0,1,0)) #(1,0)
vmat[8,] <-  c(c(0,1,0,0), c(0,0,0,1)) #(1,1)

## Ordering of cols: p(x0,y0|z1) p(x0,y1|z1) p(x1,y0|z1) ... p(x1,y1 | zp)

## {1}=1 {2}=0
vmat[9,] <-  c(c(0,0,1,0), c(1,0,0,0)) #(0,0) 
vmat[10,] <- c(c(0,0,0,1), c(1,0,0,0))#(0,1)
vmat[11,] <- c(c(0,0,1,0), c(0,1,0,0))#(1,0)
vmat[12,] <- c(c(0,0,0,1), c(0,1,0,0))#(1,1)

## Now we add two columns
##     to indicate
##     the nature of the convex hull/cone
vmat <- cbind(rep(0,4*(2^p-1)), rep(1,4*(2^p-1)),vmat)

hmat <- scdd(vmat,representation="V")

# $output
      # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 # [4,]    0    0    1    0    0    0    0    0    0     0
 # [10,]   0    0    0    1    0    0    0    0    0     0
 # [11,]   0    0    0    0    1    0    0    0    0     0
 # [12,]   0    0    0    0    0    1    0    0    0     0
 # [1,]    0    0    0    0    0    0    1    0    0     0
 # [9,]    0    0    0    0    0    0    0    1    0     0
 # [7,]    0    0    0    0    0    0    0    0    1     0
 # [8,]    0    0    0    0    0    0    0    0    0     1
 
 # [2,]    0    0    0    0    0   -1    1    1    0     1
 # [3,]    0    0    0    0   -1    0    1    1    1     0
 # [5,]    0    0    0   -1    0    0    0    1    1     1
 # [6,]    0    0   -1    0    0    0    1    0    1     1

## Note that there are only 4 additional constraints,
## not 8, as we might expect (on the basis of the conjecture
## that each observable gives rise to upper bounds)
## However, there is a simple explanation for this:
## The upper bounds on the cells in Z=1, are equivalent
## to the upper bounds on Z=0. 

