library(xtable)
library(MASS)
library(stats)
library(ggplot2)
library(gridExtra)
library(latex2exp)

source("meanFunctions.R")
source("functions.R")

reps <- 100
k <- 100
ss <- c(50,500)
n <- c()
for(i in 1:length(ss)) n <- c(n,rep(ss[i],k))

# #a,b,l1,d,l2,x,y
# B <- matrix(c(0,0,0,0,0,0,0,
#               0,0,0,0,0,0,0,
#               0,0,0,0,0,0,0,
#               1,1,1,0,0,0,0,
#               0,0,0,0,0,0,0,
#               0,1,1,0,1,0,0,
#               0,0,0,0,1,1,0),nrow=7)

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

Z <- list()
C <- list()

#a,b,l1,c,l2,x,y
#1,2,3,4,5,6,7,8
Z[[1]] <- c(1) # "(A,{B,D})"
Z[[2]] <- c(1) # "(A,{C,D})"
Z[[3]] <- c(1) # "(A,{B,C,D})"
Z[[4]] <- c(4) # "(B,C)"
Z[[5]] <- c(4) # "(B,{A,C})"
Z[[6]] <- c(4) # "(B,{C,D})"
Z[[7]] <- c(4) # "(B,{A,C,D})"
Z[[8]] <- c(5) # "(D,B)"
Z[[9]] <- c(5) # "(D,C)"
Z[[10]] <- c(5) # "(D,{B,C})"
Z[[11]] <- c(5) # "(D,{A,B})"
Z[[12]] <- c(5) # "(D,{A,C})"
Z[[13]] <- c(5) # "(D,{A,B,C})"
Z[[14]] <- c(4,5) # "({B,D},C)"
Z[[15]] <- c(4,5) # "({B,D},{A,C})"
Z[[16]] <- c(1,4) # "({A,B},C}"
Z[[17]] <- c(1,4) # "({A,B},{C,D})"
Z[[18]] <- c(1,5) # "({A,D},B)"
Z[[19]] <- c(1,5) # "({A,D},C)"
Z[[20]] <- c(1,5) # "({A,D},{B,C})"
Z[[21]] <- c(1,4,5) # "({A,B,D},C)" optimal CIS
Z[[22]] <- c() # plain OLS
Z[[23]] <- c(5) # filler to prevent bugs


#a,c,l1,b,d,l2,x,e,y
C[[1]] <- c(4,5)
C[[2]] <- c(2,5)
C[[3]] <- c(2,4,5)
C[[4]] <- c(2)
C[[5]] <- c(1,2)
C[[6]] <- c(2,5)
C[[7]] <- c(1,2,5)
C[[8]] <- c(4)
C[[9]] <- c(2)
C[[10]] <- c(2,4)
C[[11]] <- c(1,4)
C[[12]] <- c(1,2)
C[[13]] <- c(1,2,4)
C[[14]] <- c(2)
C[[15]] <- c(1,2)
C[[16]] <- c(2)
C[[17]] <- c(2,5)
C[[18]] <- c(4)
C[[19]] <- c(2)
C[[20]] <- c(2,4)
C[[21]] <- c(2)# optimal CIS
C[[22]] <- c() # plain OLS
C[[23]] <- c(2) # filler to prevent bugs

set.seed(5)
Parameters_1 <- generate_parameters(B,k) # generate Parameters
errors <- sample(c("Gaussian","uniform"),k,replace=TRUE)

k2 <- length(Z)

MSE_1 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Avar_1 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_1 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_1 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

for(j in 1:length(ss)){
  for(i in 1:k){
    Results_1 <- testGraph(n[(j-1)*k+i],reps,Parameters_1[[1]][,i],Parameters_1[[2]][,i],x,y,Z,C,B,errors[i],complex = FALSE)
    MSE_1[(j-1)*k+i,] <- Results_1[1,]
    Avar_1[(j-1)*k+i,] <- Results_1[2,]
    Var_1[(j-1)*k+i,] <- Results_1[3,]
    Bias_1[(j-1)*k+i,] <- Results_1[4,]
  }
}

aSE_1 <- sqrt(Avar_1)
RMSE_1 <- sqrt(MSE_1)

comparison_l1 <- matrix(c(21,1,
                          21,2,
                          21,3,
                          21,4,
                          21,5,
                          21,6,
                          21,7,
                          21,8,
                          21,9,
                          21,10,
                          21,11,
                          21,12,
                          21,13,
                          21,14,
                          21,15,
                          21,16,
                          21,17,
                          21,18,
                          21,19,
                          21,20
                          # 21,22
                          ),nrow=2)

names_comparisons_l1 <- c("(A,{B,D})",
                          "(A,{C,D})",
                          "(A,{B,C,D})",
                          "(B,C)",
                          "(B,{A,C})",
                          "(B,{C,D})",
                          "(B,{A,C,D})",
                          "(D,B)",
                          "(D,C)",
                          "(D,{B,C})",
                          "(D,{A,B})",
                          "(D,{A,C})",
                          "(D,{A,B,C})",
                          "({B,D},C)",
                          "{{B,D},{A,C}}",
                          "({A,B},C)",
                          "({A,B},{C,D})",
                          "({A,D},B)",
                          "({A,D},C)",
                          "({A,D},{B,C})"
                          # "OLS"
                          )

# "(B,C)" # "(B,{A,C})" # "(B,{C,D})" # "(B,{A,C,D})"
# "(D,B)" # "(D,C)" # "(D,{B,C})" # "(D,{A,B})"
# "(D,{A,C})" # "(D,{A,B,C})" # "({B,D},C)"
# "({B,D},{A,C})" # "({A,B},C}" # "({A,B},{C,D})"
# "({A,D},B)" # "({A,D},C)" # "({A,D},{B,C})"
# "({A,B,D},C)" optimal CIS
# plain OLS

comparison_l1_1 <- comparison_l1[,1:10]
names_comparisons_l1_1 <- names_comparisons_l1[1:10]
p_50_1a_l1 <- plot_MSE_comparison(RMSE_1[n==50,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1}, n=50$"),bw=.1)
p_500_1a_l1 <- plot_MSE_comparison(RMSE_1[n==500,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1}, n=500$"),bw=.1)
p_infty_1a_l1 <- plot_MSE_comparison(aSE_1[n==50,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1}"),ytitle="Asy. SD ratio",bw=.1)

pdf("example_l1_1.pdf",width=9, height=7.5)
grid.arrange(p_50_1a_l1, p_500_1a_l1, p_infty_1a_l1, ncol=1)
dev.off()


comparison_l1_2 <- comparison_l1[,11:20]
names_comparisons_l1_2 <- names_comparisons_l1[11:20]
p_50_1a_l1_2 <- plot_MSE_comparison(RMSE_1[n==50,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1}, n=50$"),bw=.1)
p_500_1a_l1_2 <- plot_MSE_comparison(RMSE_1[n==500,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1}, n=500$"),bw=.1)
p_infty_1a_l1_2 <- plot_MSE_comparison(aSE_1[n==50,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1}"),ytitle="Asy. SD ratio",bw=.1)

pdf("example_l1_2.pdf",width=9, height=7.5)
grid.arrange(p_50_1a_l1_2, p_500_1a_l1_2, p_infty_1a_l1_2, ncol=1)
dev.off()

pdf("example_l1.pdf",width=9, height=15)
grid.arrange(p_50_1a_l1, p_500_1a_l1, p_infty_1a_l1,
             p_50_1a_l1_2, p_500_1a_l1_2, p_infty_1a_l1_2, ncol=1)
dev.off()






