rm(list=ls())
### In the note below we perform the following calculation
### We derive bounds on the observable distribution under the IV model

require("rcdd")
require("gtools")
source("description_hmat.R")
p <- 4  ## dimension of Z
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

## generate a matrix represents different combinations for {X(z); z=1,...,p}
## containing 2 values
xz_combn = permutations(2,p, c(0,1), repeats.allowed = T)
## remove only all zeros and all ones
xz_combn = xz_combn[-c(1, nrow(xz_combn)),]

## generate y(x(z)=0), y(x(z)) combination
yx_combn = permutations(2, 2, c(0,1), repeats.allowed = T)
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

vmat_tmp = matrix(NA, ncol = ncol(vmat))
for(i in 1:nrow(xz_combn)){
      xz_i = xz_combn[i, ]
      tmp.mat = matrix(0, nrow = 4, ncol = ncol(vmat))
      for(j in 1: nrow(yx_combn)){
            ## go through
            ## Y(X(z)=0)=0, Y(X(z)=1)=0
            ## Y(X(z)=0)=0, Y(X(z)=1)=1
            ## Y(X(z)=0)=0, Y(X(z)=1)=1
            ## Y(X(z)=0)=1, Y(X(z)=1)=1
            yx_j = yx_combn[j, ]

            for(k in 1: length(xz_i)){
                  if(xz_i[k] == 0){ # x(z)=0
                        if(yx_j[1] == 0){ # Y(x(z)=0)=0
                              tmp.mat[j, (k-1)*4 + 1] = 1
                        }else{ # y(x(z)=0)=1
                              tmp.mat[j, (k-1)*4 + 2] = 1
                        }
                  }else{ # x(z)=1
                        if(yx_j[2] == 0){ # Y(x(z)=1)=0
                              tmp.mat[j, (k-1)*4 + 3] = 1
                        }else{ # y(x(z)=1)=1
                              tmp.mat[j, k*4] = 1
                        }
                  }
            }
      }
      vmat_tmp = rbind(vmat_tmp, tmp.mat)
}

vmat_tmp = vmat_tmp[-1,]
vmat[5:nrow(vmat), ] = vmat_tmp

# ## X(z1)=0, X(z2)=0, X(z3)=1
# vmat[5,] <-  c(c(1,0,0,0), c(1,0,0,0), c(0,0,1,0)) #(0,0)
# vmat[6,] <-  c(c(1,0,0,0), c(1,0,0,0), c(0,0,0,1)) #(0,1)
# vmat[7,] <-  c(c(0,1,0,0), c(0,1,0,0), c(0,0,1,0)) #(1,0)
# vmat[8,] <-  c(c(0,1,0,0), c(0,1,0,0), c(0,0,0,1)) #(1,1)
#
# ## Ordering of cols: p(x0,y0|z1) p(x0,y1|z1) p(x1,y0|z1) p(x1,y1|z1) ... p(x1,y1 | zp)
#
# ## X(z1)=0, X(z2)=1, X(z3)=0
# vmat[9,] <-  c(c(1,0,0,0), c(0,0,1,0), c(1,0,0,0)) #(0,0)
# vmat[10,] <- c(c(1,0,0,0), c(0,0,0,1), c(1,0,0,0))#(0,1)
# vmat[11,] <- c(c(0,1,0,0), c(0,0,1,0), c(0,1,0,0))#(1,0)
# vmat[12,] <- c(c(0,1,0,0), c(0,0,0,1), c(0,1,0,0))#(1,1)
#
# ## X(z1)=0, X(z2)=1, X(z3)=1
# vmat[13,] <- c(c(1,0,0,0), c(0,0,1,0), c(0,0,1,0)) #(0,0)
# vmat[14,] <- c(c(1,0,0,0), c(0,0,0,1), c(0,0,0,1))#(0,1)
# vmat[15,] <- c(c(0,1,0,0), c(0,0,1,0), c(0,0,1,0))#(1,0)
# vmat[16,] <- c(c(0,1,0,0), c(0,0,0,1), c(0,0,0,1))#(1,1)
#
# ## X(z1)=1, X(z2)=0, X(z3)=0


## Now we add two columns
##     to indicate
##     the nature of the convex hull/cone
vmat <- cbind(rep(0,4*(2^p-1)), rep(1,4*(2^p-1)),vmat)

hmat <- scdd(vmat,representation="V")
# View(hmat$output)
dim(hmat$output)
# 51

z3_descrip = descrip_hmat(hmat$output)

z3_descrip$eq_num
## 3
## for each z, P(x0,y0|z)+P(x0,y1|z)+P(x1,y0|z)+P(x1,y1|z)=1
z3_descrip$ineq
#      b -1  0 1 num
# [1,] 0  0 11 1   9
# [2,] 0  1  8 3   6
# [3,] 1  4  5 3   6
# [4,] 0  1  7 4   6
# [5,] 1  3  7 2  12
# [6,] 1  2 10 0   6
# [7,] 1  3  9 0   3
z3_descrip$ineq_type_index
# [[1]]
# [1]  1 10 19 28 42 44 45 46 47
#
# [[2]]
# [1]  2  4 12 16 24 27
#
# [[3]]
# [1]  3  5 11 18 20 23
#
# [[4]]
# [1]  6  9 14 15 25 26
#
# [[5]]
# [1]  7  8 13 17 21 22 29 30 34 35 38 40
#
# [[6]]
# [1] 31 32 36 37 39 41
#
# [[7]]
# [1] 33 43 48

## R1 & R7 trivial surface P(x,y|z)>=0, num = 12
## R2 & R6 : P(x1,y1|z2)+p(x1,y0|z3)<=1, num = 12
## R3 & R4 & R5: P(x0,y0|z1)+P(x1,y0|z1)+P(x1,y1|z2)+P(x0,y1|z3)+P(x1,y0|z3)<=2; (z1,z1,z2,z3,z3), (11223),(12233), num = 24
