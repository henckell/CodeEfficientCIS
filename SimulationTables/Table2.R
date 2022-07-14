library(xtable)
library(MASS)
library(stats)

source("functions.R")

testGraph <- function(epsilon,alpha,x,y,Z,W,amat,numCases){
  
  B <- amat  
  B[B!=0] <- alpha
  
  CovG <- amat_to_cov(t(B),sqrt(epsilon))	
  
  k <- length(Z)
  
  avar <- numeric(k)
  
  for(j in 1:length(Z)){
    avar[j] <- avar.twoSLS(x,y,Z[[j]],W[[j]],CovG,B[x,y])
  }
  
  return(avar)
}

set.seed(10)

k <- 6
k2 <- 3

#a,l,b,c,l2,x,y
B <- matrix(c(0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,
              0,1,1,0,0,0,0,
              0,0,0,0,0,0,0,
              1,0,0,0,1,0,0,
              0,1,0,0,1,1,0),nrow=7)

Epsilon <- matrix(runif(k*7,0.5,1),nrow=6)

Alphas <- matrix(runif(sum(B!=0)*7,0.5,2),nrow=7)
Sign2 <- matrix(rbinom(sum(B!=0)*7,1,0.5),nrow=7)
Alphas[Sign2==1] <- -Alphas[Sign2==1]

Parameters <- list(Epsilon,Alphas)

Avar <- matrix(numeric(k*k2),nrow=k)

x <- 6
y <- 7

Z <- list()
W <- list()

Z[[1]] <- c(1)
Z[[2]] <- c(1)
Z[[3]] <- c(1)

W[[1]] <- c(3,4)
W[[2]] <- c()
W[[3]] <- c(3,4)


for(i in 1:k){
  Avar[i,] <- testGraph(Epsilon[i,],Alphas[i,],x,y,Z,W,B,k)
}

print(Avar)

# testGraph(rep(1,7),rep(1,8),x,y,Z,W,B,k)

xtable(t(Avar))
