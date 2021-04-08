

blocks7.studT <-function(){
  n <- 2048
  sigma <- 7
  change.points <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-sigma*rt(n,4)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


blocks7.out <-function(){
  n <- 2048
  sigma <- 7
  change.points <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rnorm(n,mean=0,sd=sigma)
  id.out <- which(runif(n)<=.05)
  y[id.out] <- rnorm(length(id.out),mean=0,sd=4*sigma)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


blocks7.lap <-function(){
  n <- 2048
  sigma <- 4
  change.points <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rlaplace(n,s=7)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}

blocks7 <-function(){
  n <- 2048
  sigma <- 7
  change.points <- c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68, 15.37, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-c()
  for (i in meanvec){ y<- c(y,rnorm(1,i,sigma))}
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sigma=sigma))
}


###very wide teeth 


verywideteeth.25 <-function(){
  n <- 271
  sigma <- 0.25
  change.points <- c(31,61,91,121,151,181,211,241)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0,1,0,1,0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-c()
  for (i in meanvec){ y<- c(y,rnorm(1,i,sigma))}
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sigma=sigma))
}



verywideteeth.25.studT <-function(){
  n <- 271
  sigma <- 0.25
  change.points <- c(31,61,91,121,151,181,211,241)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0,1,0,1,0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-sigma*rt(n,4)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


verywideteeth.25.out <-function(){
  n <- 271
  sigma <- 0.25
  change.points <- c(31,61,91,121,151,181,211,241)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0,1,0,1,0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rnorm(n,mean=0,sd=sigma)
  id.out <- which(runif(n)<=.05)
  y[id.out] <- rnorm(length(id.out),mean=0,sd=4*sigma)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


verywideteeth.25.lap <-function(){
  n <- 271
  sigma <- 0.25
  change.points <- c(31,61,91,121,151,181,211,241)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0,1,0,1,0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rlaplace(n,s=7)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


#small wide teeth


smallwideteeth.25 <-function(){
  n <- 140
  sigma <- 0.25
  change.points <- c(31,61,91, 121)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-c()
  for (i in meanvec){ y<- c(y,rnorm(1,i,sigma))}
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sigma=sigma))
}



smallwideteeth.25.studT <-function(){
  n <- 140
  sigma <- 0.25
  change.points <- c(31,61,91, 121)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-sigma*rt(n,3)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


smallwideteeth.25.out <-function(){
  n <- 140
  sigma <- 0.25
  change.points <- c(31,61,91, 121)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rnorm(n,mean=0,sd=sigma)
  id.out <- which(runif(n)<=.1)
  y[id.out] <- rnorm(length(id.out),mean=0,sd=4*sigma)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}


smallwideteeth.25.lap <-function(){
  n <- 140
  sigma <- 0.25
  change.points <- c(31,61,91, 121)
  change.points <- c(1,change.points,n+1)
  level <- c(0, 1, 0, 1, 0)
  level <- level - level[1]
  meanvec<-rep(level,diff(change.points))
  y <-rlaplace(n,s=0.3)
  sig2hat <- var(y)
  y <- meanvec+y
  #for (i in meanvec){ y<- c(y,rnorm(1, 0.75*sigma*rt(1,5)))}
  #y <- meanvec + y
  x <- seq(1,n)/n
  sample <- cbind(x,y)
  sample <- as.data.frame(sample)
  return(list(sample=sample,change.points=change.points[-c(1,length(change.points))],level=level,sig2hat=sig2hat))
}
