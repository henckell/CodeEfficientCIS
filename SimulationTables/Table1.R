library(xtable)
library(MASS)
library(stats)

source("functions.R")

n <- 10000
reps <- 1000


testGraph <- function(epsilon,alpha,x,y,Z,W,amat,numCases){
  
  B <- amat  
  B[B!=0] <- alpha
  
  CovG <- amat_to_cov(t(B),sqrt(epsilon))	
  
  avar <- c()
  
  for(j in 1:length(Z)){
    avar[j] <- avar.twoSLS(x,y,Z[[j]],W[[j]],CovG,causalEffect(t(B),y,x,numeric(0)))
  }
  
  return(avar)
}

set.seed(20)

k <- 6
k2 <- 6

#a,c,l1,b,d,l2,x,e,y
B <- matrix(c(0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,
              0,1,0,0,0,0,0,0,0,
              1,0,1,1,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,
              0,0,1,1,0,1,0,0,0,
              0,0,0,0,0,0,1,0,0,
              0,1,0,0,0,1,1,1,0),nrow=9)

x <- 7
y <- 9

Epsilon <- matrix(runif(k*9,0.5,1),nrow=6)

Alphas <- matrix(runif(sum(B!=0)*6,0.5,2),nrow=6)
Sign2 <- matrix(rbinom(sum(B!=0)*6,1,0.5),nrow=6)
Alphas[Sign2==1] <- -Alphas[Sign2==1]

Parameters <- list(Epsilon,Alphas)
  
Avar <- matrix(numeric(k*k2),nrow=k)

Z <- list()
W <- list()

Z[[1]] <- c(5)
Z[[2]] <- c(5)
Z[[3]] <- c(5)
Z[[4]] <- c(5)
Z[[5]] <- c(1,4,5)
Z[[6]] <- c(5)

W[[1]] <- c(2)
W[[2]] <- c(1,2)
W[[3]] <- c(2,4)
W[[4]] <- c(1,2,4)
W[[5]] <- c(2)
W[[6]] <- c(4)



for(i in 1:k){
  Avar[i,] <- testGraph(Epsilon[i,],Alphas[i,],x,y,Z,W,B,k)
}

print(t(Avar))

xtable(t(Avar))
