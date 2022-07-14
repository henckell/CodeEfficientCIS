geomean <- function(data) exp(mean(log(data),na.rm=TRUE))

geomean.comp.new <- function(data, parameter, classes){
  geomean <- matrix(numeric(length(data[1,])*length(classes)),nrow=length(data[1,]))
  which.class <- list()
  for(i in 1:length(classes)) {
    which.class[[i]] <- which(parameter==classes[i])
    geomean[,i] <- apply(data[which.class[[i]],],2,geomean) 
  }
  colnames(geomean) <- classes
  rownames(geomean) <- c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG","posspaX est. CPDAG","paX est. CPDAG"
                         ,"O est. CPDAG","empty with 0")
  return(geomean)
}


mean.comp <- function(data, parameter, classes,names=TRUE){
  mean <- matrix(numeric(length(data[1,])*length(classes)),nrow=length(data[1,]))
  which.class <- list()
  for(i in 1:length(classes)) {
    which.class[[i]] <- which(parameter==classes[i])
    mean[,i] <- apply(data[which.class[[i]],],2,mean,na.rm=TRUE) 
  }
  colnames(mean) <- classes
  if(names==TRUE){
    rownames(mean) <- c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG"
                        ,"posspaX est. CPDAG","paX est. CPDAG","O est. CPDAG","empty with 0")
    }
  return(mean)
}


median.comp <- function(data, parameter, classes,names=TRUE){
  median <- matrix(numeric(length(data[1,])*length(classes)),nrow=length(data[1,]))
  which.class <- list()
  for(i in 1:length(classes)) {
    which.class[[i]] <- which(parameter==classes[i])
    median[,i] <- apply(data[which.class[[i]],],2,median,na.rm=TRUE) 
  }
  colnames(median) <- classes
  if(names==TRUE){
    rownames(median) <- c("empty","paX DAG","O DAG","paX est. CPDAG","O est. CPDAG","empty with 0")
  }
  return(median)
}

geomean.plot <- function(data,parameter,classes,names = classes,name,
                         lim=c(0.0001,10),pos="left",ylab="geomean MSE"){
  geomean <- geomean.comp.new(data,parameter,classes)
  
  plot(y=t(geomean),log="y", x=rep(1:length(classes),length(geomean[,1])), ylim= lim, xlim = c(1,length(classes))
       ,xaxt ="n",xlab=name, ylab=ylab)  
  axis(1,at=1:length(classes),labels = names)
  for(i in 1:length(geomean[,1])){
    lines(geomean[i,],col = i,lty=i)
  }
  if(pos=="left"){
    legend("topleft",legend=c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG","posspaX est. CPDAG","paX est. CPDAG"
                              ,"O est. CPDAG","empty with 0"),col=1:length(geomean[,1]),pch=1,lty=1:length(geomean[,1]))
  }
  if(pos=="right"){
    legend("topright",legend=c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG","posspaX est. CPDAG","paX est. CPDAG"
                               ,"O est. CPDAG","empty with 0"),col=1:length(geomean[,1]),pch=1,lty=1:length(geomean[,1]))
  }
}

median.plot <- function(data,parameter,classes,names = classes,name
                        ,lim=c(0.0001,10),pos="left",title="",name2="median MSE"){
  median <- median.comp(data,parameter,classes)
  
  plot(y=t(median),log="y",x=rep(1:length(classes),length(median[,1]))
       , ylim= lim, xlim = c(1,length(classes)), 
       xaxt ="n",xlab=name, ylab=name2,main=title)  
  axis(1,at=1:length(classes),labels = names)
  for(i in 1:length(median[,1])){
    lines(median[i,],col = i,lty=i)
  }
  if(pos=="left"){
    legend("topleft",legend=c("empty","paX DAG","O DAG","paX est. CPDAG","O est. CPDAG","empty with 0")
           ,col=1:length(median[,1]),pch=1,lty=1:length(median[,1]))
  }
  if(pos=="right"){
    legend("topright",legend=c("empty","paX DAG","O DAG","paX est. CPDAG","O est. CPDAG","empty with 0")
           ,col=1:length(median[,1]),pch=1,lty=1:length(median[,1]))
  }
}

mean.plot <- function(data,parameter,classes,names = classes,name,lim=c(0.0001,10000),pos="left",names2=TRUE,title="",
                      name2="mean MSE"){
  mean <- mean.comp(data,parameter,classes,names2)
  
  plot(y=t(mean),log="y",x=rep(1:length(classes),length(mean[,1])), ylim= lim, xlim = c(1,length(classes)), 
       xaxt ="n",xlab=name, ylab=name2,main=title)  
  axis(1,at=1:length(classes),labels = names)
  for(i in 1:length(mean[,1])){
    lines(mean[i,],col = i,lty=i)
  }
  if(pos=="left"){
    legend("topleft",legend=c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG","posspaX est. CPDAG","paX est. CPDAG"
                              ,"O est. CPDAG","empty with 0"),col=1:length(mean[,1]),pch=1,lty=1:length(mean[,1]))
  }
  if(pos=="right"){
    legend("topright",legend=c("empty","paX DAG","O DAG","posspaX true CPDAG","paX true CPDAG","O true CPDAG","posspaX est. CPDAG","paX est. CPDAG"
                               ,"O est. CPDAG","empty with 0"),col=1:length(mean[,1]),pch=1,lty=1:length(mean[,1]))
  }
}

