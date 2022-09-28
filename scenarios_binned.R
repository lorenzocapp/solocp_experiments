blocks7.out.mod.r <-function(){
  n <- 2048/2
  sigma <- 7
  change.points <- c(204,  472, 820, 1332, 1658)/2
  change.points <- c(0,change.points,n)
  change.points.here <- change.points/n
  level <- c(0, 14.64,  -7.32, 3.29, 19.03, 0)
  level <- level - level[1]
  y <- rep(NA,n)
  x <- rep(NA,n)
  for (i in 1:n){
    u<-runif(1)
    x[i] <- u
    loc.mean<- min(which(u<=change.points.here))-1
    if (runif(1)>.1){
        y[i] <-level[loc.mean]+rnorm(1,mean=0,sd=sigma)
     } else {
        y[i] <-level[loc.mean]+rnorm(1,mean=0,sd=4*sigma)
     }
  }
  sig2hat <- var(y)
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


blocks7.studT.r <-function(){
  n <- 2048/2
  sigma <- 7
  change.points <- c(204,  472, 820, 1332, 1658)/2
  change.points <- c(0,change.points,n)
  change.points.here <- change.points/n
  level <- c(0, 14.64,  -7.32, 3.29, 19.03, 0)
  level <- level - level[1]
  y <-sigma*rt(n,4)
  x <- rep(NA,n)
  for (i in 1:n){
    u<-runif(1)
    x[i] <- u
    loc.mean<- min(which(u<=change.points.here))-1
    y[i] <-level[loc.mean]+y[i]
  }
  sig2hat <- var(y)
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


blocks7.lap.r <-function(){
  n <- 2048/2
  sigma <- 7
  change.points <- c(204,  472, 820, 1332, 1658)/2
  change.points <- c(0,change.points,n)
  change.points.here <- change.points/n
  level <- c(0, 14.64,  -7.32, 3.29, 19.03, 0)
  level <- level - level[1]
  y <-rlaplace(n,s=9)
  x <- rep(NA,n)
  for (i in 1:n){
    u<-runif(1)
    x[i] <- u
    loc.mean<- min(which(u<=change.points.here))-1
    y[i] <-level[loc.mean]+y[i]
  }
  sig2hat <- var(y)
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


blocks7.r <-function(){
  n <- 2048/2
  sigma <- 7
  change.points <- c(204,  472, 820, 1332, 1658)/2
  change.points <- c(0,change.points,n)
  change.points.here <- change.points/n
  level <- c(0, 14.64,  -7.32, 3.29, 19.03, 0)
  level <- level - level[1]
  y <-rnorm(n,0,sigma)
  x <- rep(NA,n)
  for (i in 1:n){
    u<-runif(1)
    x[i] <- u
    loc.mean<- min(which(u<=change.points.here))-1
    y[i] <-level[loc.mean]+y[i]
  }
  sig2hat <- var(y)
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


