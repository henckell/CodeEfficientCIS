library(xtable)
library(MASS)
library(stats)
library(ggplot2)
library(gridExtra)
library(latex2exp)

source("meanFunctions.R")
source("functions.R")

reps <- 100
k <- 1000
ss <- c(20,500)
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

set.seed(60)
Parameters_2 <- generate_parameters(B,k) # generate Parameters
errors <- sample(c("Gaussian","uniform"),k,replace=TRUE)

k2 <- length(Z)

MSE_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Avar_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

for(j in 1:length(ss)){
  for(i in 1:k){
    Results_2 <- testGraph(n[(j-1)*k+i],reps,Parameters_2[[1]][,i],Parameters_2[[2]][,i],x,y,Z,C,B,errors[i])
    MSE_2[(j-1)*k+i,] <- Results_2[1,]
    Avar_2[(j-1)*k+i,] <- Results_2[2,]
    Var_2[(j-1)*k+i,] <- Results_2[3,]
    Bias_2[(j-1)*k+i,] <- Results_2[4,]
  }
}

aSE_2 <- sqrt(Avar_2)
RMSE_2 <- sqrt(MSE_2)

comparison_l1 <- matrix(c(18,1,
                       18,2,
                       18,3,
                       18,4,
                       18,5,
                       18,6,
                       18,7,
                       18,8,
                       18,9,
                       18,10,
                       18,11,
                       # 18,12,
                       18,13,
                       18,14,
                       18,15,
                       18,16,
                       18,17,
                       18,19),nrow=2)

names_comparisons_l1 <- c("{B,C}",
                       "{B,{A,C}}",
                       "{B,{C,D}}",
                       "{B,{A,C,D}}",
                       "{D,B}", 
                       "{D,C}",
                       "{D,{B,C}}",
                       "{D,{A,B}}",
                       "{D,{A,C}}",
                       "{D,{A,B,C}}",
                       "{{B,D},C}",
                       # "{{B,D},{A,C}}",
                       "{{A,B},C}",
                       "{{A,B},{C,D}}",
                       "{{A,D},B}",
                       "{{A,D},C}",
                       "{{A,D},{B,C}}",
                       "OLS")

# "(B,C)" # "(B,{A,C})" # "(B,{C,D})" # "(B,{A,C,D})"
# "(D,B)" # "(D,C)" # "(D,{B,C})" # "(D,{A,B})"
# "(D,{A,C})" # "(D,{A,B,C})" # "({B,D},C)" 
# "({B,D},{A,C})" # "({A,B},C}" # "({A,B},{C,D})"
# "({A,D},B)" # "({A,D},C)" # "({A,D},{B,C})"
# "({A,B,D},C)" optimal CIS
# plain OLS

comparison_l1_1 <- comparison_l1[,1:8]
names_comparisons_l1_1 <- names_comparisons_l1[1:8]
p_20_1a_l1 <- plot_MSE_comparison(RMSE_2[n==20,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1a}, n=20$"))
p_500_1a_l1 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1a}, n=500$"))
p_infty_1a_l1 <- plot_MSE_comparison(aSE_2[n==20,],comparison_l1_1,names_comparisons_l1_1,title=TeX("$G_{1a}"),ytitle="Asy. SD ratio")

pdf("plots/example_l1_1.pdf",width=9, height=7.5)
grid.arrange(p_20_1a_l1, p_500_1a_l1, p_infty_1a_l1, ncol=1)
dev.off()


comparison_l1_2 <- comparison_l1[,9:17]
names_comparisons_l1_2 <- names_comparisons_l1[9:17]
p_20_1a_l1_2 <- plot_MSE_comparison(RMSE_2[n==20,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1a}, n=20$"))
p_500_1a_l1_2 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1a}, n=500$"))
p_infty_1a_l1_2 <- plot_MSE_comparison(aSE_2[n==20,],comparison_l1_2,names_comparisons_l1_2,title=TeX("$G_{1a}"),ytitle="Asy. SD ratio")

pdf("plots/example_l1_2.pdf",width=9, height=7.5)
grid.arrange(p_20_1a_l1_2, p_500_1a_l1_2, p_infty_1a_l1_2, ncol=1)
dev.off()

pdf("plots/example_l1.pdf",width=9, height=15)
grid.arrange(p_20_1a_l1, p_500_1a_l1, p_infty_1a_l1,
             p_20_1a_l1_2, p_500_1a_l1_2, p_infty_1a_l1_2, ncol=1)
dev.off()





   
