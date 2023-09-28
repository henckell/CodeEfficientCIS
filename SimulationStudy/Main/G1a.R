library(xtable)
library(MASS)
library(stats)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(SteinIV)
library(ivmodel)


setwd("~/GitHub/efficientCIS")
source("HelperFunctions/functions.R")

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

# Example 2 CIS:
#a,b,l1,c,l2,x,y
#1,2,3,4,5,6,7,8

Z[[1]] <- c(4) # "(B,C)"
Z[[2]] <- c(4) # "(B,{A,C})"
Z[[3]] <- c(4) # "(B,{C,D})"
Z[[4]] <- c(4) # "(B,{A,C,D})"
Z[[5]] <- c(5) # "(D,B)"
Z[[6]] <- c(5) # "(D,C)"
Z[[7]] <- c(5) # "(D,{B,C})"
Z[[8]] <- c(5) # "(D,{A,B})"
Z[[9]] <- c(5) # "(D,{A,C})"
Z[[10]] <- c(5) # "(D,{A,B,C})"
Z[[11]] <- c(4,5) # "({B,D},C)"
Z[[12]] <- c(4,5) # "({B,D},{A,C})"
Z[[13]] <- c(1,4) # "({A,B},C}"
Z[[14]] <- c(1,4) # "({A,B},{C,D})"
Z[[15]] <- c(1,5) # "({A,D},B)"
Z[[16]] <- c(1,5) # "({A,D},C)"
Z[[17]] <- c(1,5) # "({A,D},{B,C})"
Z[[18]] <- c(1,4,5) # "({A,B,D},C)" optimal CIS
Z[[19]] <- c() # plain OLS
Z[[20]] <- c(5) # filler to prevent bugs

C[[1]] <- c(2)
C[[2]] <- c(1,2)
C[[3]] <- c(2,5)
C[[4]] <- c(1,2,5)
C[[5]] <- c(4)
C[[6]] <- c(2)
C[[7]] <- c(2,4)
C[[8]] <- c(1,4)
C[[9]] <- c(1,2)
C[[10]] <- c(1,2,4)
C[[11]] <- c(2)
C[[12]] <- c(1,2)
C[[13]] <- c(2)
C[[14]] <- c(2,5)
C[[15]] <- c(4)
C[[16]] <- c(2)
C[[17]] <- c(2,4)
C[[18]] <- c(2)# optimal CIS
C[[19]] <- c() # plain OLS
C[[20]] <- c(2) # filler to prevent bugs

set.seed(5)
Parameters <- generate_parameters(B,k) # generate Parameters
errors <- sample(c("Gaussian","uniform"),k,replace=TRUE)

k2 <- length(Z)

MSE_tsls <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Avar_tsls <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_tsls <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_tsls <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

MSE_liml <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_liml <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_liml <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

MSE_jive <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_jive <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_jive <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

for(j in 1:length(ss)){
  for(i in 1:k){
    Results <- testGraph(n[(j-1)*k+i],reps,Parameters[[1]][,i],Parameters[[2]][,i],x,y,Z,C,B,errors[i])
    MSE_tsls[(j-1)*k+i,] <- Results$tsls[1,]
    Avar_tsls[(j-1)*k+i,] <- Results$tsls[2,]
    Var_tsls[(j-1)*k+i,] <- Results$tsls[3,]
    Bias_tsls[(j-1)*k+i,] <- Results$tsls[4,]

    MSE_liml[(j-1)*k+i,] <- Results$liml[1,]
    Var_liml[(j-1)*k+i,] <- Results$liml[2,]
    Bias_liml[(j-1)*k+i,] <- Results$liml[3,]

    MSE_jive[(j-1)*k+i,] <- Results$jive[1,]
    Var_jive[(j-1)*k+i,] <- Results$jive[2,]
    Bias_jive[(j-1)*k+i,] <- Results$jive[3,]
  }
}

aSE_tsls_1 <- sqrt(Avar_tsls)
RMSE_tsls_1 <- sqrt(MSE_tsls)
RMSE_jive_1 <- sqrt(MSE_jive)
RMSE_liml_1 <- sqrt(MSE_liml)

comparison_1 <- matrix(c(
                       18,1,
                       # 18,2,
                       # 18,3,
                       # 18,4,
                       18,5,
                       18,6,
                       18,7,
                       18,10,
                       18,11,
                       # 15,14,
                       18,78,
                       18,58,
                       18,19
                       ),nrow=2)

names_comparisons_1 <- c(
                       "(B,C)",
                       # "(B,{A,C})",
                       # "(B,{C,D})",
                       # "(B,{A,C,D})",
                       "(D,B)",
                       "(D,C)",
                       "(D,{B,C})",
                       "(D,{A,B,C})",
                       "({B,D},C)",
                       # "{{B,C},{A}}",
                       "LIML",
                       "JIVE",
                       "OLS")

RMSE_joint_1 <- cbind(RMSE_tsls_1,aSE_tsls_1,RMSE_jive_1,RMSE_liml_1)

p_50_1a <- plot_MSE_comparison(RMSE_joint_1[n==50,],comparison_1,names_comparisons_1,title=TeX("$G_{1}, n=50$"),bw=0.1)
p_500_1a <- plot_MSE_comparison(RMSE_joint_1[n==500,],comparison_1,names_comparisons_1,title=TeX("$G_{1}, n=500$"),bw=0.1)
p_infty_1a <- plot_MSE_comparison(aSE_tsls_1[n==50,],comparison_1,names_comparisons_1,title=c(""),ytitle="Asy. SD ratio")

pdf("plots/example2.pdf",width=9, height=7.5)
grid.arrange(p_50_1a, p_500_1a, p_infty_1a, ncol=1)
dev.off()

write.table(RMSE_joint_1, file = "Example1_MSE_100.txt", sep = ",", quote = FALSE, row.names = F)

comparison_1_2 <- matrix(c(18,19,
                         18,58,
                         18,78),nrow=2)

names_comparisons_1_2 <- c("OLS","JIVE","LIML")


p_50_1a_2 <- plot_MSE_comparison(RMSE_joint_1[n==50,],comparison_1_2,names_comparisons_1_2,title=TeX("$G_{1}, n=50$"))
p_500_1a_2 <- plot_MSE_comparison(RMSE_joint_1[n==500,],comparison_1_2,names_comparisons_1_2,title=TeX("$G_{1}, n=500$"),ytitle="",bw=0.1)

pdf("example1and2.pdf",width=9, height=7.5)
grid.arrange(p_50_1a, p_500_1a, p_50_2a, p_500_2a, ncol=1)
dev.off()

pdf("method_comp.pdf",width=9, height=3.75)
grid.arrange(p_50_1a_2, p_500_1a_2, p_50_2a_2, p_500_2a_2, ncol=2, nrow=2)
dev.off()

# full names
# c("({B},em)",
#   "({B},{A})",
#   "({B},{C})",
#   "({B},{A,C})",
#   "({C},em)",  # 5
#   "({C},{A})",
#   "({C},{B})",
#   "({C},{A,B})",
#   "({A,B},em)",
#   "({A,B},{C})",  # 10
#   "({A,C},em)",
#   "({A,C},{B})",
#   "({B,C},em)",
#   "({B,C},{A})",
#   "({A,B,C},em)",  # 15
#   "em")
