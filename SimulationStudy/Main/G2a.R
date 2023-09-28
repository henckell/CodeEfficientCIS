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



#a,b,c,l,d,x,y
B <- matrix(c(0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,
              0,1,0,1,1,0,0,
              0,0,1,1,0,1,0),nrow=7)

x <- 6
y <- 7

Z <- list()
C <- list()

#a,b,c,l,d,x,y
#1,2,3,4,5,6,7
Z[[1]] <- c(1)
Z[[2]] <- c(1)
Z[[3]] <- c(1)
Z[[4]] <- c(1)
Z[[5]] <- c(2)
Z[[6]] <- c(2)
Z[[7]] <- c(2)
Z[[8]] <- c(2)
Z[[9]] <- c(2)
Z[[10]] <- c(2)
Z[[11]] <- c(2)
Z[[12]] <- c(2)
Z[[13]] <- c(5) # "({D},em)"
Z[[14]] <- c(5) # "({D},\{A\})"
Z[[15]] <- c(5) # "({D},\{B\})"
Z[[16]] <- c(5) # "({D},\{C\})"
Z[[17]] <- c(5) # "({D},\{A,B\})"
Z[[18]] <- c(5) # "({D},\{A,C\})"
Z[[19]] <- c(5) # "({D},\{B,C\})"
Z[[20]] <- c(5) # "({D},\{A,B,C\})"
Z[[21]] <- c(1,2)
Z[[22]] <- c(1,2)
Z[[23]] <- c(1,2)
Z[[24]] <- c(1,2)
Z[[25]] <- c(1,5)
Z[[26]] <- c(1,5)
Z[[27]] <- c(1,5)
Z[[28]] <- c(1,5)
Z[[29]] <- c(2,5) # optimal instrument
Z[[30]] <- c(2,5)
Z[[31]] <- c(2,5) # optimal CIS
Z[[32]] <- c(2,5)
Z[[33]] <- c(1,2,5) # optimal instrument
Z[[34]] <- c(1,2,5) # optimal CIS 2
Z[[35]] <- c() # plain OLS
Z[[36]] <- c(1) # filler to prevent bugs

C[[1]] <- c() # single instrument A
C[[2]] <- c(2)
C[[3]] <- c(3)
C[[4]] <- c(2,3)
C[[5]] <- c() # single instrument B
C[[6]] <- c(1)
C[[7]] <- c(3)
C[[8]] <- c(5)
C[[9]] <- c(1,3)
C[[10]] <- c(1,5)
C[[11]] <- c(3,5)
C[[12]] <- c(1,3,5)
C[[13]] <- c() # single instrument D
C[[14]] <- c(1)
C[[15]] <- c(2)
C[[16]] <- c(3) # plain OLS
C[[17]] <- c(1,2) # filler to prevent bugs
C[[18]] <- c(1,3)
C[[19]] <- c(2,3)
C[[20]] <- c(1,2,3)
C[[21]] <- c() # instrument \{A,B\}
C[[22]] <- c(3)
C[[23]] <- c(5)
C[[24]] <- c(3,5)
C[[25]] <- c() # instrument \{A,D\}
C[[26]] <- c(2)
C[[27]] <- c(3)
C[[28]] <- c(2,3)
C[[29]] <- c() # optimal instrument \{B,D\}
C[[30]] <- c(1)
C[[31]] <- c(3) # optimal CIS
C[[32]] <- c(1,3)
C[[33]] <- c() # instrument \{A,B,D\}
C[[34]] <- c(3) # optimal CIS 2
C[[35]] <- c() # plain OLS
C[[36]] <- c(2) # filler to prevent bugs

set.seed(60)
Parameters <- generate_parameters(B,k) # generate Parameters
errors <- sample(c("Gaussian","uniform"),k,replace=TRUE)
# c("Gaussian","uniform")

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

aSE_tsls_2 <- sqrt(Avar_tsls)
RMSE_tsls_2 <- sqrt(MSE_tsls)
RMSE_jive_2 <- sqrt(MSE_jive)
RMSE_liml_2 <- sqrt(MSE_liml)
RMSE_joint_2 <- cbind(RMSE_tsls_2,aSE_tsls_2,RMSE_jive_2,RMSE_liml_2)

comparison_2 <- matrix(c(31,13,
                       31,14,
                       31,15,
                       31,16,
                       # 31,18,
                       31,20,
                       # 31,21,
                       # 31,25,
                       31,29,
                       # 31,33,
                       31,139,
                       31,103,
                       31,35
                       ),nrow=2)

names_comparisons_2 <- c("(D,\u2205)",
                       "(D,A)",
                       "(D,B)",
                       "(D,C)",
                       # "(D,{A,C})",
                       "(D,{A,B,C})",
                       # "({A,B},\u2205)",
                       # "({A,D},\u2205)",
                       "({B,D},\u2205)",
                       # "{{A,B,D},em}",
                       "LIML",
                       "JIVE",
                       "OLS"
                       )

p_50_2a <- plot_MSE_comparison(RMSE_joint_2[n==50,],comparison_2,names_comparisons_2,title=TeX("$G_{2}, n=50$"),bw=.1)
p_500_2a <- plot_MSE_comparison(RMSE_joint_2[n==500,],comparison_2,names_comparisons_2,title=TeX("$G_{2}, n=500$"),bw=.1)

pdf("example2.pdf",width=9, height=7.5)
grid.arrange(p_50_2a, p_500_2a, p_infty_2a, ncol=1)
dev.off()

grid.arrange(p_50_2a, p_500_2a, p_infty_2a, ncol=1)

write.table(RMSE_joint_2, file = "Example2_MSE_100.txt", sep = ",", quote = FALSE, row.names = F)

comparison_2_2 <- matrix(c(31,35,
                           31,103,
                           31,139),nrow=2)

names_comparisons_2_2 <- c("OLS","JIVE","LIML")

p_50_2a_2 <- plot_MSE_comparison(RMSE_joint_2[n==50,],comparison_2_2,names_comparisons_2_2,title=TeX("$G_{2}, n=50$"))
p_500_2a_2 <- plot_MSE_comparison(RMSE_joint_2[n==500,],comparison_2_2,names_comparisons_2_2,title=TeX("$G_{2}, n=500$"),ytitle="",bw=.1)



comparison_alt_2 <- matrix(c(31,13,
                         31,14,
                         31,15,
                         31,16,
                         31,18,
                         31,20,
                         31,21,
                         31,25,
                         31,29
                         # 31,33,
),nrow=2)

names_alt_2 <- c("(D,\u2205)",
                         "(D,A)",
                         "(D,B)",
                         "(D,C)",
                         "(D,{A,C})",
                         "(D,{A,B,C})",
                         "({A,B},\u2205)",
                         "({A,D},\u2205)",
                         "({B,D},\u2205)"
                         # "{{A,B,D},em}",
)


p_50_2a_jive <- plot_MSE_comparison(RMSE_jive_2[n==50,],comparison_alt_2,names_alt_2,title=TeX("$G_{2}, n=50, JIVE$"),bw=0.1)
p_500_2a_jive <- plot_MSE_comparison(RMSE_jive_2[n==500,],comparison_alt_2,names_alt_2,title=TeX("$G_{2}, n=500, JIVE$"),bw=0.1)
p_50_2a_liml <- plot_MSE_comparison(RMSE_liml_2[n==50,],comparison_alt_2,names_alt_2,title=TeX("$G_{2}, n=50, LIML$"),bw=0.1)
p_500_2a_liml <- plot_MSE_comparison(RMSE_liml_2[n==500,],comparison_alt_2,names_alt_2,title=TeX("$G_{2}, n=500, LIML$"),bw=0.1)

pdf("limljive_2.pdf",width=9, height=3.75)
grid.arrange(p_50_2a_liml, p_500_2a_liml,p_50_2a_jive, p_500_2a_jive, nrow=4)
dev.off()
