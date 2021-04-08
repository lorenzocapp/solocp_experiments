ebpiece <- function(y, sig2=1, alpha, v, lambda, M) {
  
  n <- length(y)
  burn <- round(0.5 * M)
  MM <- M + burn
  BB <- matrix(0, nrow=MM, ncol=n - 1)
  theta <- matrix(0, nrow=MM, ncol=n)
  B <- sample(0:1, n - 1, replace=TRUE)
  b <- sum(B) + 1
  bb <- numeric(MM)
  BB[1,] <- B
  bb[1] <- b
  ii <- integer(MM)
  ii[1] <- i <- iuB <- nuB <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(B, y)
  vv <- 1 + alpha * v / sig2 
  sq.vv <- sqrt(vv)
  lprior <- -lambda * (1:n - 1) * log(n) - lchoose(n - 1, 1:n - 1)
  lpost <- lprior[b] - alpha * out[[1]]$sse / 2 / sig2 - b * log(vv) / 2
  for(m in 1:MM) {
    
    B.new <- rprop(B)
    b.new <- sum(B.new) + 1
    i.new <- compare.blocks(as.matrix(BB[iuB,]), B.new)
    if(i.new == 0) o.new <- get.lm.stuff(B.new, y) else o.new <- out[[i.new]]
    lpost.new <- lprior[b.new] - alpha * o.new$sse / 2 / sig2 - b.new * log(vv) / 2
    if(runif(1) <= exp(lpost.new - lpost)) {
      
      B <- B.new
      b <- b.new
      lpost <- lpost.new
      if(i.new == 0) {
        
        nuB <- nuB + 1
        i <- nuB
        iuB <- c(iuB, m)
        out[[nuB]] <- o.new
        
      } else i <- i.new
      
    }
    BB[m,] <- B
    bb[m] <- b
    ii[m] <- i
    theta[m,] <- rep(out[[i]]$theta.B + sqrt(v / vv / out[[i]]$B.split) * rnorm(b), out[[i]]$B.split)
    
  }
  return(list(theta=theta[-(1:burn),], B=BB[-(1:burn),], b=bb[-(1:burn)]))
  
  
}



ebpiece_new <- function(y, sig2=1, alpha, v, lambda, M,B) {
  
  n <- length(y)
  burn <- round(0.5 * M)
  MM <- M + burn
  BB <- matrix(0, nrow=MM, ncol=n - 1)
  theta <- matrix(0, nrow=MM, ncol=n)
  #B <- sample(0:1, n - 1, replace=TRUE)
  
  b <- sum(B) + 1
  bb <- numeric(MM)
  BB[1,] <- B
  bb[1] <- b
  ii <- integer(MM)
  ii[1] <- i <- iuB <- nuB <- 1
  out <- list()
  out[[1]] <- get.lm.stuff(B, y)
  vv <- 1 + alpha * v / sig2 
  sq.vv <- sqrt(vv)
  lprior <- -lambda * (1:n - 1) * log(n) - lchoose(n - 1, 1:n - 1)
  lpost <- lprior[b] - alpha * out[[1]]$sse / 2 / sig2 - b * log(vv) / 2
  for(m in 1:MM) {
    
    B.new <- rprop(B)
    b.new <- sum(B.new) + 1
    i.new <- compare.blocks(as.matrix(BB[iuB,]), B.new)
    if(i.new == 0) o.new <- get.lm.stuff(B.new, y) else o.new <- out[[i.new]]
    lpost.new <- lprior[b.new] - alpha * o.new$sse / 2 / sig2 - b.new * log(vv) / 2
    if(runif(1) <= exp(lpost.new - lpost)) {
      
      B <- B.new
      b <- b.new
      lpost <- lpost.new
      if(i.new == 0) {
        
        nuB <- nuB + 1
        i <- nuB
        iuB <- c(iuB, m)
        out[[nuB]] <- o.new
        
      } else i <- i.new
      
    }
    BB[m,] <- B
    bb[m] <- b
    ii[m] <- i
    theta[m,] <- rep(out[[i]]$theta.B + sqrt(v / vv / out[[i]]$B.split) * rnorm(b), out[[i]]$B.split)
    
  }
  return(list(theta=theta[-(1:burn),], B=BB[-(1:burn),], b=bb[-(1:burn)]))
  
  
}

# Auxiliary functions

rprop <- function(B) {
  
  i <- sample(seq_along(B), 1)
  B[i] <- (B[i] + 1) %% 2
  return(B)
  
}


compare.blocks <- function(B.old, B.new) {
  
  h <- function(v) sum(abs(v - B.new))
  o <- apply(B.old, 1, h)
  if(all(o > 0)) return(0) else return(which(o == 0)[1])
  
}


get.lm.stuff <- function(B, y) {
  
  if(all(B == 0)) {
    
    theta.B <- mean(y)
    sse <- sum((y - theta.B)**2)
    B.split <- length(y)
    
  }
  else {
    
    x <- as.factor(c(1, 1 + cumsum(B)))
    o <- lm(y ~ x - 1)
    sse <- sum(o$residuals**2)
    theta.B <- o$coefficients
    B.split <- as.numeric(table(x))
    
  }
  return(list(sse=sse, theta.B=theta.B, B.split=B.split))
  
}


# Do examples

ebpiece.example <- function(theta, sig2, v, lambda, M) {
  
  n <- length(theta)
  y <- theta + sqrt(sig2) * rnorm(n)
  o <- ebpiece(y, sig2, 0.99, v, lambda, M)
  plot(y, cex=0.5, col="gray")
  lines(theta)
  points(apply(o$theta, 2, mean), pch=19, cex=0.5)
  tb <- table(o$b)
  plot(tb / sum(tb), xlab="|B|", ylab="Probability")
  return(o)
  
}


dist_change_points =  function(Shat,S0)
{
  
  if(length(Shat)==0)
  {return(Inf)}
  
  if(length(S0)==0)
  {return(-Inf)}
  
  temp =rep(0,length(S0))
  for(j in 1:length(S0))
  {
    temp[j] = min(abs(S0[j] - Shat))
  }
  return( max(temp) )
}

N.vec_gen <- function(est.cp,true.cp){
  if (length(est.cp)==-0){
    return(N.vec <- c(length(true.cp),rep(0,6)))
  }
  
  N.vec <- rep(0,7)
  for (j in true.cp) {
    i.min<-which(abs(est.cp-j)==min(abs(est.cp-j)))[1]
    diff <- est.cp[i.min]-j
    if      (diff<=-3) {N.vec[1]=N.vec[1]+1
    } else if (diff==-2) {N.vec[2]=N.vec[2]+1
    } else if (diff==-1) {N.vec[3]=N.vec[3]+1
    } else if (diff==0) {N.vec[4]=N.vec[4]+1
    } else if (diff==1) {N.vec[5]=N.vec[5]+1
    } else if (diff==2) {N.vec[6]=N.vec[6]+1
    } else {N.vec[7]=N.vec[7]+1}
  }
  return(N.vec/length(true.cp)) 
} 

N.vec_gen.inv <- function(est.cp,true.cp){
  if (length(est.cp)==0){
    return(N.vec <- c(length(true.cp),rep(0,6)))
  }
  
  N.vec <- rep(0,7)
  for (j in est.cp) {
    i.min<-which(abs(true.cp-j)==min(abs(true.cp-j)))[1]
    diff <- true.cp[i.min]-j
    if      (diff<=-3) {N.vec[1]=N.vec[1]+1
    } else if (diff==-2) {N.vec[2]=N.vec[2]+1
    } else if (diff==-1) {N.vec[3]=N.vec[3]+1
    } else if (diff==0) {N.vec[4]=N.vec[4]+1
    } else if (diff==1) {N.vec[5]=N.vec[5]+1
    } else if (diff==2) {N.vec[6]=N.vec[6]+1
    } else {N.vec[7]=N.vec[7]+1}
  }
  return(N.vec/length(est.cp)) 
} 


basad <- function(X,Y,Z0,B0 = NULL,returnB = 0, K,del,sig,nburn = 1000, 
                  niter = 5000,nsplit = 1){
  p = dim(X)[2]-1; n = dim(X)[1]
  choicep <-function(x){
    return(x - K + qnorm(0.9)*sqrt(x*(1- x/p)))
  }
  cp = uniroot(choicep, c(1,K))$root
  
  #s1  =  sig*0.1*p^{2 + del}/n;
  s1  =  max(1, 0.01*p^{2 + del}/n);
  s0 = (0.1*sig)/(n);      ##Informative prior params
  
  
  ### nsplit is the number of iterative conditional normals
  vsize = (p+1) %/%nsplit 
  G = t(X)%*%X
  gib.beta <- function(X, Y, B, sig, Z){
    T1 = Z*s1 + (1-Z)*s0
    vec = seq(1:vsize)
    Xtmp = X
    for(s in 1:nsplit){
      svec = (s-1)*vsize +vec
      COV = (G[svec,svec] + diag(1/T1[svec]))
      ec = eigen(COV)
      COVsq = ec$vectors %*% diag(1/sqrt(ec$values)) %*% t(ec$vectors)
      B[svec] = COVsq %*% (COVsq %*% (t(X[,svec]) %*% Y - G[svec, -svec] %*% B[-svec]) + sqrt(sig)*rnorm(vsize))
    }
    if( p+1 > nsplit*vsize){
      svec = (nsplit*vsize +1):(p+1)
      COV = (G[svec,svec] + diag(1/T1[svec]))
      ec = eigen(COV)
      COVsq = ec$vectors %*% diag(1/sqrt(ec$values)) %*% t(ec$vectors)
      B[svec] = COVsq %*% (COVsq %*% (t(X[,svec]) %*% Y - G[svec, -svec] %*% B[-svec]) + sqrt(sig)*rnorm(length(svec)))
    }
    return(B)
  }
  
  p1 = 10^-4; p2 = 10^-4                 ## Paramerers of IG distribution for Error variance
  
  gib.sig <- function(X, B, Z){
    T1 = Z*s1 + (1-Z)*s0
    res = 1/rgamma(1, p1+ (n*0.5)+(p*0.5), p2+(0.5* t(Y - X %*% B)%*% (Y - X%*%B))+(0.5*t(B)%*%diag(1/T1)%*%B))
    return(res)
  }
  
  gib.z <- function(B, sig, pr){
    A = array(0, p+1)
    s = seq(1:(p+1))
    prob = sapply(s, function(j) pr* dnorm(B[j], 0, sqrt(sig*s1))/ (pr* dnorm(B[j], 0, sqrt(sig*s1))+ (1-pr)* dnorm(B[j], 0, sqrt(sig*s0))))
    tmp = (runif((p+1)) - prob)
    Z = (tmp <0);
    
    if(sum(Z) > n/2) {
      indz = which(Z ==1)
      Z[indz[which(prob[Z] < rev(sort(prob[Z]))[round(n/2)])]] = 0}
    return(Z)
  }
  
  ## Gibbs algorithm
  
  B = array(0, c(niter,p+1))
  Z = array(0, c(niter,p+1))
  Z[1,] = Z0
  if(is.null(B0) ==FALSE) {B[1,] = B0}
  err.var = array(0, niter)
  pr = cp/p
  err.var[1] = 1
  
  
  for(brn in 2:nburn){
    
    B[brn,]  = gib.beta(X,Y,  B[(brn-1),], err.var[(brn-1)], Z[(brn-1),])
    err.var[brn] = gib.sig(X, B[brn,], Z[(brn-1),])
    Z[brn,] = gib.z(B[brn,], err.var[brn], pr)
    Z[brn,1] = 1
    #if(brn %%1000 ==0) print(brn)
  }
  
  
  B[1,] = B[nburn,]
  Z[1,] = Z[nburn,]
  pr = cp/p
  err.var[1] = err.var[nburn]
  BB = B[1,]
  ZZ = ZZ1 = apply(Z[1:nburn,],2,mean)
  
  for(itr in 2:niter){
    B[itr,]  = gib.beta(X,Y, B[(itr-1),], err.var[(itr-1)], Z[(itr-1),])
    err.var[itr] = gib.sig(X, B[itr,], Z[(itr-1),])
    Z[itr,] = gib.z(B[itr,], err.var[itr], pr)
    Z[itr,1] = 1
    #if(itr %%1000 ==0) print(itr)
  }
  if(returnB==0) {return(Z)}
  if(returnB==1) {return(list(Z=Z,B=B))}
}

