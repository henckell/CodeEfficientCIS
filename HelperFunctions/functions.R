amat_to_cov <- function(amat, eps_sd=NULL, eps_var=1) {
  "Derive the covariance matrix from a given adjacency matrix and the error standard deviances.
  Equation by Bollen (pg. 80)."
  
  # get number of intended variables
  n_var = ncol(amat)
  # error means and variances are either given as vectors or scalars. create vector if needed
  if (!is.null(eps_sd)) {eps_var = eps_sd^2}
  if (length(eps_var) == 1) {eps_var = rep(eps_var[[1]], n_var)}
  
  I = diag(n_var)
  ci = amat - I
  ci_inv = solve(ci)
  rownames(ci_inv) = rownames(amat)
  colnames(ci_inv) = colnames(amat)
  cmat = ci_inv %*% diag(eps_var) %*% t(ci_inv)
}

conCov <- function(A, B, C=NULL,cov){
  "Calculates the conditional covariance from a given covariance matrix and two variable names.
  If C is an empty set (NULL), then it just returns the normal covariance."
  
  # trivial input with empty conditioning set
  if (is.null(C)) {return(cov[A,B])}
  
  (cov[A,B] - cov[A,C] %*% solve(cov[C,C]) %*% cov[C,B])
}

avar.twoSLS <- function(x,y,z,c,Cov,tau){
  # compute asymptotic variance according to our equation
  
  if(is.null(z)&&is.null(c)){
    Y.X <-  conCov(y,y,x,Cov)
    return(Y.X/Cov[x,x])
  }
  
  X.C <- conCov(x,x,c,Cov)
  X.ZC <- conCov(x,x,union(z,c),Cov)
  
  hatY.C <- conCov(y,y,c,Cov) - 2 * tau * conCov(x,y,c,Cov) + tau^2 * X.C
  
  aVarTau <-   hatY.C/(X.C- X.ZC)
  
  return(aVarTau)
}

twoSLS <- function(x,y,z,c=NULL,dat){
  if(is.null(z)&&is.null(c)){
    return(solve(t(dat[,x])%*%dat[,x])%*%(t(dat[,x])%*%dat[,y]))
  }
  
  if(is.null(z)){
    a <- union(x,c)
    return((solve(t(dat[,a])%*%dat[,a])%*%(t(dat[,a])%*%dat[,y]))[1])
  }
  
  a <- union(x,c)
  b <- union(z,c)
  
  gamma <- solve(t(dat[,b]) %*% dat[,b]) %*% (t(dat[,b]) %*% dat[,a])
  
  (solve(t(gamma) %*% t(dat[,b]) %*% dat[,b] %*% gamma) %*% (t(gamma) %*% t(dat[,b]) %*% dat[,y]))[1]
}


generate_parameters <- function(B,k){
  
  num_nodes <- nrow(B)
  Epsilon <- matrix(runif(k*num_nodes,0.1,1),nrow=num_nodes)
  Sign <- matrix(rbinom(k*num_nodes,1,0.5),nrow=num_nodes)
  
  num_edges <- sum(B==1)
  Alphas <- matrix(runif(k*num_edges,0.1,2),nrow=num_edges)
  Sign2 <- matrix(rbinom(k*num_edges,1,0.5),nrow=num_edges)
  Alphas[Sign2==1] <- -Alphas[Sign2==1]
  
  Parameters <- list(Epsilon=Epsilon,Alphas=Alphas)
  return(Parameters)
}

causalEffect <- function (B, y, x, z) {
  # requires topolically ordered B
  p <- ncol(B)
  vec <- matrix(0, p, 1)
  vec[x] <- 1
  if (y - x > 1) {
    for (i in (x + 1):y){
      if(sum(i==z)==0){vec[i] <- B[i, ] %*% vec} 
      else{vec[i]=0}
    }
    return(vec[y])
  }
  else {
    return(B[y, x])
  }
}

## asssumes that forbidden projection has already been applied to amat

testGraph <- function(n,reps,epsilon,alpha,x,y,Z,C,amat,error="Gaussian"){
  
  B <- amat  
  B[B!=0] <- alpha
  
  CovG <- amat_to_cov(t(B),sqrt(epsilon))	
  
  k <- length(Z)
  avar <- numeric(k)
  
  betas <- matrix(numeric(reps*k),ncol=k)
  
  
  for(i in 1:reps){
    if(error=="Gaussian"){
    Data <- mvrnorm(n,rep(0,length(epsilon)),CovG)}
    else{
      if(error=="uniform"){
        Data <- DataGenUnif(n,t(B),epsilon)
      }
      if(error=="student"){
        Data <- DataGenUnif(n,t(B),epsilon)
      }
    }
    for(j in 1:length(Z)){
      betas[i,j] <- twoSLS(x,y,Z[[j]],C[[j]],Data)
    }
  }
  
  bias <- apply(abs(betas-causalEffect(t(B),y,x,numeric(0))),2,mean)
  var <- apply(betas,2,stats::var)
  MSE <- apply((betas-causalEffect(t(B),y,x,numeric(0)))^2,2,mean)
  
  for(j in 1:length(Z)){
    avar[j] <- avar.twoSLS(x,y,Z[[j]],C[[j]],CovG,causalEffect(t(B),y,x,numeric(0)))
  }
  
  return(rbind(MSE,avar,var,bias,n))
}

plot_MSE_comparison <- function(MSE, comparison, name_comparison,title=c(""),xtitle=c(""),ytitle=c("RMSE ratio")){
  
  x <- y <- gm <- med <- c()
  for(i in 1:ncol(comparison)){
    y <- c(y,MSE[,comparison[1,i]]/MSE[,comparison[2,i]])
    x <- c(x,rep(i,nrow(MSE)))
    gm <- cbind(gm,geomean(c(MSE[,comparison[1,i]]/MSE[,comparison[2,i]])))
    med <- cbind(med,median(c(MSE[,comparison[1,i]]/MSE[,comparison[2,i]]),na.rm = TRUE))
  }
  
  geo.df <- data.frame(y=t(gm),x=1:ncol(comparison))
  med.df <- data.frame(y=t(med),x=1:ncol(comparison))
  
  df <- as.data.frame(x=factor(x),y=y)
  
  p <-  ggplot(df, aes(x=factor(x),y=y),fill=factor(x)) + geom_violin() +
    # geom_abline(slope=0,intercept=1,color = "blue") + 
    ylim(0,2) +
    ggtitle(title) +
    ylab(ytitle) + xlab(xtitle) +
    scale_x_discrete(labels=name_comparison) +
    geom_point(shape=15,data=geo.df,aes(x=x,y=y), colour = 'black',inherit.aes=FALSE) 
    # geom_point(shape=15,data=med.df,aes(x=x,y=y), colour = 'black',inherit.aes=FALSE)
  
  p
}

## uniform errors
DataGenUnif<-function(sample.size,B,var){
  limit <- sqrt(3*var)
  p <- length(B[1,])
  X<-matrix(numeric(sample.size*p),ncol=p)
  X[,1] = runif(sample.size,-limit[1],limit[1])
  X[,2] = B[2,1] %*% X[,1] + runif(sample.size,-limit[2],limit[2])
  if(p>=3){
    for(i in 3:p){
      X[,i] <-X[,1:(i-1)]  %*% t(t(B[i,1:(i-1)])) + runif(sample.size,-limit[i],limit[i])
    }
  }
  return(X)
}

## student-t distributed errors
DataGenStu<-function(sample.size,B,df=10,var){
  if(length(df)==1) df <- rep(df,nrow(B))
  p <- length(B[1,])
  X<-matrix(numeric(sample.size*p),ncol=p)
  X[,1] = rt(sample.size,df[1]) * sqrt(var[1] * (df[1]-2)/df[1])
  X[,2] = B[2,1] %*% X[,1] + rt(sample.size,df[2]) * sqrt(var[2] * (df[2]-2)/df[2])
  if(p>=3){
    for(i in 3:p){
      X[,i] <-X[,1:(i-1)]  %*% t(t(B[i,1:(i-1)])) + rt(sample.size,df[i]) * sqrt(var[i] * (df[i]-2)/df[i])
    }
  }
  return(X)
}
