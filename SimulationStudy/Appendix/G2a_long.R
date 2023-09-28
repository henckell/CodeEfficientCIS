library(xtable)
library(MASS)
library(stats)
library(ggplot2)
library(gridExtra)
library(latex2exp)

#### take the square root of the MSE
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
C[[16]] <- c(3)
C[[17]] <- c(1,2)
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
Parameters_1 <- generate_parameters(B,k) # generate Parameters
errors <- sample(c("Gaussian","uniform"),k,replace=TRUE)
# c("Gaussian","uniform")

k2 <- length(Z)

MSE_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Avar_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Var_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)
Bias_2 <- matrix(numeric(length(ss)*k*k2),nrow=length(ss)*k)

for(j in 1:length(ss)){
  for(i in 1:k){
    Results_2 <- testGraph(n[(j-1)*k+i],reps,Parameters_1[[1]][,i],Parameters_1[[2]][,i],x,y,Z,C,B
                           ,error=errors[i],complex = FALSE)
    MSE_2[(j-1)*k+i,] <- Results_2[1,]
    Avar_2[(j-1)*k+i,] <- Results_2[2,]
    Var_2[(j-1)*k+i,] <- Results_2[3,]
    Bias_2[(j-1)*k+i,] <- Results_2[4,]
  }
}

aSE_2 <- sqrt(Avar_2)
RMSE_2 <- sqrt(MSE_2)

comparison_l2 <- rbind(c(rep(31,34)),c(1:30,32:35))


names_comparisons_l2 <- c("(A,\u2205)",
                          "(A,B)",
                          "(A,C)",
                          "(A,{B,C})",
                          "(B,\u2205)",
                          "(B,A)",
                          "(B,C)",
                          "(B,D)",
                          "(B,{A,C})",
                          "(B,{A,D})",
                          "(B,{C,D})",
                          "(B,{A,B,D})",
                          "(D,\u2205)",
                          "(D,A)",
                          "(D,B)",
                          "(D,C)",
                          "(D,{A,B})",
                          "(D,{A,C})",
                          "(D,{B,C})",
                          "(D,{A,B,C})",
                          "({A,B},\u2205)",
                          "({A,B},C)",
                          "({A,B},D)",
                          "({A,B},{C,D})",
                          "({A,D},\u2205)",
                          "({A,D},B)",
                          "({A,D},C)",
                          "({A,D},{B,C})",
                          "({B,D},\u2205)",
                          "({B,D},A)",
                          "({B,D},{A,C})",
                          "({A,B,D},\u2205)",
                          "({A,B,D},C)",
                          "OLS"
                          )

comparison_l2_1 <- comparison_l2[,1:8]
names_comparisons_l2_1 <- names_comparisons_l2[1:8]
p_50_2a_l2_1 <- plot_MSE_comparison(RMSE_2[n==50,],comparison_l2_1,names_comparisons_l2_1,title=TeX("$G_{2}, n=50$"),bw=0.1)
p_500_2a_l2_1 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l2_1,names_comparisons_l2_1,title=TeX("$G_{2}, n=500$"),bw=0.1)
p_infty_2a_l2_1 <- plot_MSE_comparison(aSE_2[n==50,],comparison_l2_1,names_comparisons_l2_1,title=TeX("$G_{2}$"),ytitle="Asy. SD ratio",bw=0.1)

pdf("example_l2_1.pdf",width=9, height=7.5)
grid.arrange(p_50_2a_l2_1, p_500_2a_l2_1, p_infty_2a_l2_1, ncol=1)
dev.off()

comparison_l2_2 <- comparison_l2[,9:16]
names_comparisons_l2_2 <- names_comparisons_l2[9:16]
p_50_2a_l2_2 <- plot_MSE_comparison(RMSE_2[n==50,],comparison_l2_2,names_comparisons_l2_2,title=TeX("$G_{2}, n=50$"),bw=0.1)
p_500_2a_l2_2 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l2_2,names_comparisons_l2_2,title=TeX("$G_{2}, n=500$"),bw=0.1)
p_infty_2a_l2_2 <- plot_MSE_comparison(aSE_2[n==50,],comparison_l2_2,names_comparisons_l2_2,title=TeX("$G_{2}$"),ytitle="Asy. SD ratio",bw=0.1)

pdf("example_l2_2.pdf",width=9, height=7.5)
grid.arrange(p_50_2a_l2_2, p_500_2a_l2_2, p_infty_2a_l2_2, ncol=1)
dev.off()

comparison_l2_3 <- comparison_l2[,17:25]
names_comparisons_l2_3 <- names_comparisons_l2[17:25]
p_50_2a_l2_3 <- plot_MSE_comparison(RMSE_2[n==50,],comparison_l2_3,names_comparisons_l2_3,title=TeX("$G_{2}, n=50$"),bw=0.1)
p_500_2a_l2_3 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l2_3,names_comparisons_l2_3,title=TeX("$G_{2}, n=500$"),bw=0.1)
p_infty_2a_l2_3 <- plot_MSE_comparison(aSE_2[n==50,],comparison_l2_3,names_comparisons_l2_3,title=TeX("$G_{2}$"),ytitle="Asy. SD ratio",bw=0.1)

pdf("example_l2_3.pdf",width=9, height=7.5)
grid.arrange(p_50_2a_l2_3, p_500_2a_l2_3, p_infty_2a_l2_3, ncol=1)
dev.off()

comparison_l2_4 <- comparison_l2[,c(25:33)]
names_comparisons_l2_4 <- names_comparisons_l2[c(25:33)]
p_50_2a_l2_4 <- plot_MSE_comparison(RMSE_2[n==50,],comparison_l2_4,names_comparisons_l2_4,title=TeX("$G_{2}, n=50$"),bw=0.1)
p_500_2a_l2_4 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l2_4,names_comparisons_l2_4,title=TeX("$G_{2}, n=500$"),bw=0.1)
p_infty_2a_l2_4 <- plot_MSE_comparison(aSE_2[n==50,],comparison_l2_4,names_comparisons_l2_4,title=TeX("$G_{2}$"),ytitle="Asy. SD ratio",bw=0.1)

pdf("example_l2_4.pdf",width=9, height=7.5)
grid.arrange(p_50_2a_l2_4, p_500_2a_l2_4, p_infty_2a_l2_4, ncol=1)
dev.off()


# comparison_l2_5 <- as.matrix(comparison_l2[,33])
# names_comparisons_l2_5 <- names_comparisons_l2[33]
# p_50_2a_l2_5 <- plot_MSE_comparison(RMSE_2[n==50,],comparison_l2_5,names_comparisons_l2_5,title=TeX("$G_{2a}, n=20$"))
# p_500_2a_l2_5 <- plot_MSE_comparison(RMSE_2[n==500,],comparison_l2_5,names_comparisons_l2_5,title=TeX("$G_{2a}, n=500$"))
# p_infty_2a_l2_5 <- plot_MSE_comparison(aSE_2[n==50,],comparison_l2_5,names_comparisons_l2_5,title=TeX("$G_{2a}$"),ytitle="Asy. SD ratio")
#
# pdf("example_l2_5.pdf",width=9, height=7.5)
# grid.arrange(p_50_2a_l2_5, p_500_2a_l2_5, p_infty_2a_l2_5, ncol=1)
# dev.off()


pdf("example_l2_a.pdf",width=9, height=15)
grid.arrange(p_50_2a_l2_1, p_500_2a_l2_1, p_infty_2a_l2_1,
             p_50_2a_l2_2, p_500_2a_l2_2, p_infty_2a_l2_2, ncol=1)
dev.off()

pdf("example_l2_b.pdf",width=9, height=15)
grid.arrange(p_50_2a_l2_3, p_500_2a_l2_3, p_infty_2a_l2_3,
             p_50_2a_l2_4, p_500_2a_l2_4, p_infty_2a_l2_4, ncol=1)
dev.off()
