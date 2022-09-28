#######################################
###### FUNCTIONS TREND FILTERING ######
#######################################

solo_spike_slab_single<-function(sample,tau2,tau2.spike,tau2.slab,q,sigma2){
  n<-length(sample$x)
  grid <- sort(sample$x)
  n.grid <- length(grid)-1
  n1<-rep(1,n.grid+1)
  y<-sample[order(sample$x),2]
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau2/(tau2+sigma2) #this is going down
  lower.local.means<-y[n.grid+1] #this is going down
  for (j1 in (n-1):1){
    new.lower<-tau2*((n-j1+1)-sum(lower))^2/(tau2*((n-j1+1)-(sum(lower)))+sigma2)
    new.lower.mean<-(sum(y[j1:n])-sum(lower.local.means*lower))/((n-j1+1)-sum(lower))
    lower<-c(new.lower,lower)
    lower.local.means<-c(new.lower.mean,lower.local.means)
  }
  
  sum.lower<-rev(cumsum(rev(lower)))
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-rev(cumsum(rev(mean.disc)))
  marg.mean <- c(sum.mean.disc[2:length(sum.mean.disc)],0)
  
  
  marg.weight <- c(sum.lower[2:length(sum.lower)],0)
  #Define M and GAMMA and Y matrices
  sum.inv<-seq(n,1)
  
  GAM<-matrix(NA,nrow=n,ncol=n)
  M<-matrix(NA,nrow=n,ncol=n)
  Y<-matrix(NA,nrow=n,ncol=n)
  
  sum.y<-rev(cumsum(rev(y)))
  
  #param
  GAM[,1] <- 1
  Y[,1] <- sum.y[1]-marg.mean
  M[,1]<-1/((sum.inv[1] - marg.weight)*GAM[,1]+sigma2*tau2^(-1))
  GAM[2:n,2]<-1-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]
  Y[2:n,2]<-sum.y[2]-marg.mean[2:n]-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]*GAM[2:n,1]*Y[2:n,1]
  #weight
  w.spike<-sqrt(tau2.spike^(-1)/(lower[1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.spike^(-1)))#maybe wrong but it does not matter
  w.slab<-sqrt(tau2.slab^(-1)/(lower[1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.slab^(-1)))#maybe wrong but it does not matter
  w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/((sum.inv[2] - marg.weight[2])*GAM[2,1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*Y[2,2]^2/((sum.inv[2] - marg.weight[2])*GAM[2,1]+tau2.spike^(-1))))
  w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/((sum.inv[2] - marg.weight[2])*GAM[2,1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*Y[2,2]^2/((sum.inv[2] - marg.weight[2])*GAM[2,1]+tau2.slab^(-1))))
  for (i in 3:(n-1)){
    #param
    M[(i-1):n,(i-1)]<-1/((sum.inv[(i-1)] - marg.weight[(i-1):n])*GAM[(i-1):n,(i-1)]+sigma2*tau2^(-1))
    GAM[i:n,i]<-1-(sum.inv[i] - marg.weight[i:n])*rowSums(M[i:n,1:(i-1)]*GAM[i:n,1:(i-1)]^2)
    Y[i:n,i]<-sum.y[i]-marg.mean[i:n]-(sum.inv[i] - marg.weight[i:n])*rowSums(M[i:n,1:(i-1)]*GAM[i:n,1:(i-1)]*Y[i:n,1:(i-1)])
    #weight
    w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/((sum.inv[i] - marg.weight[i])*GAM[i,i-1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*Y[i,i]^2/((sum.inv[i] - marg.weight[i])*GAM[i,i-1]+tau2.spike^(-1))))
    w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/((sum.inv[i] - marg.weight[i])*GAM[i,i-1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*Y[i,i]^2/((sum.inv[i] - marg.weight[i])*GAM[i,i-1]+tau2.slab^(-1))))
    
  }
  #param
  M[(n-1):n,(n-1)]<-1/((sum.inv[(n-1)] - marg.weight[(n-1):n])*GAM[(n-1):n,(n-1)]+sigma2*tau2^(-1))
  GAM[n,n]<-1-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]^2)
  Y[n,n]<-sum.y[i]-marg.mean[n]-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]*Y[n,1:(n-1)])
  #weight
  w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/((sum.inv[n] - marg.weight[n])*GAM[n-1,n-1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*Y[n,n]^2/((sum.inv[n] - marg.weight[n])*GAM[n-1,n-1]+tau2.spike^(-1))))
  w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/((sum.inv[n] - marg.weight[n])*GAM[n-1,n-1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*Y[n,n]^2/((sum.inv[n] - marg.weight[n])*GAM[n-1,n-1]+tau2.slab^(-1))))
  
  
  
  
  ratio<-q*w.slab/(q*w.slab+(1-q)*w.spike)
  id.inf <- which(w.slab==Inf)
  ratio[id.inf] <- 1
  
  
  return(list(w.slab=w.slab,
              w.spike=w.spike,
              ratio=ratio))
  
  
}





subset_changepoints <- function(ratio,del=5){
  first <- which(ratio >=.5)
  n.c <- length(first)
  
  low <- first - del
  #id.low <- which(low[2:(n.c)]<= first[1:(n.c-1)])
  #id.up <- which(up[1:(n.c-1)]>= first[2:(n.c)]) you do not need both sides: it is symmetric!
  id.split <-which((low[2:(n.c)]<= first[1:(n.c-1)])==FALSE)
  
  i.low.bl <-first[c(1,id.split+1)]
  i.up.bl <- first[c(id.split,n.c)]
  
  
  #sweep through the block to pick a change point.
  change.points <- c()
  for (i in 1: length(i.low.bl)){
    f <- which(ratio[i.low.bl[i]:i.up.bl[i]]==max(ratio[i.low.bl[i]:i.up.bl[i]]))
    change.points <- c(change.points, i.low.bl[i]+ f[ceiling(length(f)/2)] -1)
  }
  #remove change.points right at the beginning and at the end
  change.points<-change.points[change.points>5 & change.points<(length(ratio)-5)]
  return(change.points)
}

solo_spike_slab_grid_free<-function(sample,tau2,tau2.spike,tau2.slab,q,sigma2){
  n01<-table(sample$x)
  grid <- sort(sample$x)
  n.grid <- length(grid)-1
  idn<-which(as.character(grid) %in% names(n01))#this finds which elements are not observed
  n1<-rep(0,n.grid+1)
  n1[idn]<-n01
  yy<-tapply(sample$y,sample$x, FUN=sum)
  group.sums.y<-rep(0,n.grid+1)
  group.sums.y[idn]<-yy
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau2*n1[n.grid+1]^2/(n1[n.grid+1]*tau2+sigma2) #this is going down
  for (j1 in (n.grid):1){
    new.lower<-tau2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))^2/(tau2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))+sigma2)
    lower<-c(new.lower,lower)
  }
  sum.lower<-rev(cumsum(rev(lower)))
  marg.weigths<-matrix(rep(c(sum.lower[3:length(sum.lower)],0),n.grid+1),ncol=n.grid+1)
  
  
  #Find the weighted means
  s<-data.frame(sample)
  local<-aggregate(s[,2],list(s$x),mean)$x
  local.means<-rep(0,n.grid+1) #CHECK IF IT IS OK TO PUT 0
  local.means[idn]<-local
  
  #Compute group local means
  lower.local.means<-local.means[n.grid+1] #this is going down
  for (j1 in (n.grid):1){
    new.lower.mean<-(sum(local.means[j1:(n.grid+1)]*n1[j1:(n.grid+1)])-sum(lower.local.means*lower[(j1+1):(n.grid+1)]))/(sum(n1[j1:(n.grid+1)])-(sum(lower[(j1+1):(n.grid+1)])))
    lower.local.means<-c(new.lower.mean,lower.local.means)
  }
  
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-rev(cumsum(rev(mean.disc)))
  marg.means<-matrix(rep(c(sum.mean.disc[3:length(sum.mean.disc)],0),n.grid+1),ncol=n.grid+1)
  
  #Define M and GAMMA and Y matrices
  sum.inv<-rev(cumsum(rev(n1)))
  sum.matrix<-matrix(rep(sum.inv,n.grid),ncol=length(sum.inv),byrow=TRUE)
  GAM<-matrix(1,nrow=n.grid,ncol=n.grid+1)
  M<-matrix(0,nrow=n.grid,ncol=n.grid+1)
  sum.y<-rev(cumsum(rev(group.sums.y)))
  Y<-matrix(rep(sum.y,n.grid),ncol=length(sum.y),byrow=TRUE)
  
  #corrected n sums
  sum.matrix.prime<-sum.matrix-marg.weigths
  Y.prime<-Y-marg.means
  
  M[,1]<-tau2/(tau2*sum.matrix.prime[,1]*GAM[,1]+sigma2)
  GAM[,2]<-GAM[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]^2
  Y.prime[,2]<-Y.prime[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]*Y.prime[,1]
  for (i in 3:(n.grid+1)){
    M[,(i-1)]<-tau2/(tau2*sum.matrix.prime[,(i-1)]*GAM[,(i-1)]+sigma2)
    GAM[,i]<-GAM[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]^2)
    Y.prime[,i]<-Y.prime[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]*Y.prime[,1:(i-1)])
  }
  
  
  #Find posterior means
  id<-seq(2,n.grid+1)
  post.means<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2^(-1))
  for (i in id){
    post.means<-c(post.means,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2^(-1)))
  }
  
  id<-seq(2,n.grid+1)
  post.means.spike<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2.spike^(-1))
  post.means.slab<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2.slab^(-1))
  w.spike<-sqrt(tau2.spike^(-1)/(lower[1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.spike^(-1)))
  w.slab<-sqrt(tau2.slab^(-1)/(lower[1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.slab^(-1)))
  for (i in id){
    post.means.spike<-c(post.means.spike,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1)))
    post.means.slab<-c(post.means.slab,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1)))
    w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*Y.prime[i-1,i]^2/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1))))
    w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*Y.prime[i-1,i]^2/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1))))
  }
  
  
  ratio<-q*w.slab/(q*w.slab+(1-q)*w.spike)
  id.inf <- which(w.slab==Inf)
  ratio[id.inf] <- 1
  
  post.means.spike.and.slab<-post.means.spike
  post.means.spike.and.slab[ratio>=0.5]<-post.means.slab[ratio>=0.5]
  
  return(list(post.means=post.means,
              post.means.slab=post.means.slab,
              post.means.spike=post.means.spike,
              post.means.spike.and.slab=post.means.spike.and.slab,
              w.slab=w.slab,
              w.spike=w.spike,
              ratio=ratio))
  
  
}



solo_spike_slab<-function(sample,tau2,tau2.spike,tau2.slab,q,sigma2,grid){
  n01<-table(sample$x)
  idn<-which(as.character(grid) %in% names(n01))#this finds which elements are not observed
  n1<-rep(0,n.grid+1)
  n1[idn]<-n01
  yy<-tapply(sample$y,sample$x, FUN=sum)
  group.sums.y<-rep(0,n.grid+1)
  group.sums.y[idn]<-yy
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau2*n1[n.grid+1]^2/(n1[n.grid+1]*tau2+sigma2) #this is going down
  for (j1 in (n.grid):1){
    new.lower<-tau2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))^2/(tau2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))+sigma2)
    lower<-c(new.lower,lower)
  }
  sum.lower<-rev(cumsum(rev(lower)))
  marg.weigths<-matrix(rep(c(sum.lower[3:length(sum.lower)],0),n.grid+1),ncol=n.grid+1)
  
  
  #Find the weighted means
  s<-data.frame(sample)
  local<-aggregate(s[,2],list(s$x),mean)$x
  local.means<-rep(0,n.grid+1) #CHECK IF IT IS OK TO PUT 0
  local.means[idn]<-local
  
  #Compute group local means
  lower.local.means<-local.means[n.grid+1] #this is going down
  for (j1 in (n.grid):1){
    new.lower.mean<-(sum(local.means[j1:(n.grid+1)]*n1[j1:(n.grid+1)])-sum(lower.local.means*lower[(j1+1):(n.grid+1)]))/(sum(n1[j1:(n.grid+1)])-(sum(lower[(j1+1):(n.grid+1)])))
    lower.local.means<-c(new.lower.mean,lower.local.means)
  }
  
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-rev(cumsum(rev(mean.disc)))
  marg.means<-matrix(rep(c(sum.mean.disc[3:length(sum.mean.disc)],0),n.grid+1),ncol=n.grid+1)
  
  #Define M and GAMMA and Y matrices
  sum.inv<-rev(cumsum(rev(n1)))
  sum.matrix<-matrix(rep(sum.inv,n.grid),ncol=length(sum.inv),byrow=TRUE)
  GAM<-matrix(1,nrow=n.grid,ncol=n.grid+1)
  M<-matrix(0,nrow=n.grid,ncol=n.grid+1)
  sum.y<-rev(cumsum(rev(group.sums.y)))
  Y<-matrix(rep(sum.y,n.grid),ncol=length(sum.y),byrow=TRUE)
  
  #corrected n sums
  sum.matrix.prime<-sum.matrix-marg.weigths
  Y.prime<-Y-marg.means
  
  M[,1]<-tau2/(tau2*sum.matrix.prime[,1]*GAM[,1]+sigma2)
  GAM[,2]<-GAM[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]^2
  Y.prime[,2]<-Y.prime[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]*Y.prime[,1]
  for (i in 3:(n.grid+1)){
    M[,(i-1)]<-tau2/(tau2*sum.matrix.prime[,(i-1)]*GAM[,(i-1)]+sigma2)
    GAM[,i]<-GAM[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]^2)
    Y.prime[,i]<-Y.prime[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]*Y.prime[,1:(i-1)])
  }
  
  
  #Find posterior means
  id<-seq(2,n.grid+1)
  post.means<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2^(-1))
  for (i in id){
    post.means<-c(post.means,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2^(-1)))
  }
  
  id<-seq(2,n.grid+1)
  post.means.spike<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2.spike^(-1))
  post.means.slab<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2.slab^(-1))
  w.spike<-sqrt(tau2.spike^(-1)/(lower[1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.spike^(-1)))
  w.slab<-sqrt(tau2.slab^(-1)/(lower[1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.slab^(-1)))
  for (i in id){
    post.means.spike<-c(post.means.spike,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1)))
    post.means.slab<-c(post.means.slab,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1)))
    w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*Y.prime[i-1,i]^2/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.spike^(-1))))
    w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*Y.prime[i-1,i]^2/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2.slab^(-1))))
  }
  
  
  ratio<-q*w.slab/(q*w.slab+(1-q)*w.spike)
  
  
  post.means.spike.and.slab<-post.means.spike
  post.means.spike.and.slab[ratio>=0.5]<-post.means.slab[ratio>=0.5]
  
  return(list(post.means=post.means,
              post.means.slab=post.means.slab,
              post.means.spike=post.means.spike,
              post.means.spike.and.slab=post.means.spike.and.slab,
              w.slab=w.slab,
              w.spike=w.spike,
              ratio=ratio))
  
  
}





conj_tf_change_points<-function(sample,tau,ratio,grid){
  
  old.grid<-grid
  new.sample<-sample
  
  id<-unique(c(1,which(ratio>0.5),n.grid+1))
  for (i in 1:(length(id)-1)){
    new.sample$x[new.sample$x>=grid[id[i]] & new.sample$x<(grid[id[i+1]])]<-grid[id[i]]
  }
  new.sample$x[new.sample$x==grid[id[i+1]]]<-grid[id[i]]
  
  sample<-new.sample
  
  grid<-sort(unique(sample$x))
  n.grid<-length(grid)-1
  
  n01<-table(sample$x)
  idn<-which(as.character(grid) %in% names(n01))#this finds which elements are not observed
  n1<-rep(0,n.grid+1)
  n1[idn]<-n01
  yy<-tapply(sample$y,sample$x, FUN=sum)
  group.sums.y<-rep(0,n.grid+1)
  group.sums.y[idn]<-yy
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau^2*n1[n.grid+1]^2/(n1[n.grid+1]*tau^2+1) #this is going down
  for (j1 in (n.grid):1){
    new.lower<-tau^2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))^2/(tau^2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))+1)
    lower<-c(new.lower,lower)
  }
  sum.lower<-rev(cumsum(rev(lower)))
  marg.weigths<-matrix(rep(c(sum.lower[3:length(sum.lower)],0),n.grid+1),ncol=n.grid+1)
  
  
  #Find the weighted means
  s<-data.frame(sample)
  local<-aggregate(s[,2],list(s$x),mean)$x
  local.means<-rep(0,n.grid+1) #CHECK IF IT IS OK TO PUT 0
  local.means[idn]<-local
  
  #Compute group local means
  lower.local.means<-local.means[n.grid+1] #this is going down
  for (j1 in (n.grid):1){
    new.lower.mean<-(sum(local.means[j1:(n.grid+1)]*n1[j1:(n.grid+1)])-sum(lower.local.means*lower[(j1+1):(n.grid+1)]))/(sum(n1[j1:(n.grid+1)])-(sum(lower[(j1+1):(n.grid+1)])))
    lower.local.means<-c(new.lower.mean,lower.local.means)
  }
  
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-rev(cumsum(rev(mean.disc)))
  marg.means<-matrix(rep(c(sum.mean.disc[3:length(sum.mean.disc)],0),n.grid+1),ncol=n.grid+1)
  
  #Define M and GAMMA and Y matrices
  sum.inv<-rev(cumsum(rev(n1)))
  sum.matrix<-matrix(rep(sum.inv,n.grid),ncol=length(sum.inv),byrow=TRUE)
  GAM<-matrix(1,nrow=n.grid,ncol=n.grid+1)
  M<-matrix(0,nrow=n.grid,ncol=n.grid+1)
  sum.y<-rev(cumsum(rev(group.sums.y)))
  Y<-matrix(rep(sum.y,n.grid),ncol=length(sum.y),byrow=TRUE)
  
  #corrected n sums
  sum.matrix.prime<-sum.matrix-marg.weigths
  Y.prime<-Y-marg.means
  
  M[,1]<-tau^2/(tau^2*sum.matrix.prime[,1]*GAM[,1]+1)
  GAM[,2]<-GAM[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]^2
  Y.prime[,2]<-Y.prime[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]*Y.prime[,1]
  for (i in 3:(n.grid+1)){
    M[,(i-1)]<-tau^2/(tau^2*sum.matrix.prime[,(i-1)]*GAM[,(i-1)]+1)
    GAM[,i]<-GAM[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]^2)
    Y.prime[,i]<-Y.prime[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]*Y.prime[,1:(i-1)])
  }
  
  
  #Find posterior means
  id<-seq(2,n.grid+1)
  post.means<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2^(-1))
  for (i in id){
    post.means<-c(post.means,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2^(-1)))
  }
  
  
  #Fix the error
  #post.means[which(ratio<=0.5)]<-0
  which(old.grid  %in% grid)
  corr.post.means<-rep(0,length(old.grid))
  corr.post.means[which(old.grid  %in% grid)]<-post.means
  return(post.means=corr.post.means)
  
  
}




#to fix: try to remove the part in which I am putting zeros at the end: the way we 
#handle unobserved grid points seems wrong also in the above  
# this is the correction. If you remove it you see the problme
#I am talking about post.means[which(tf$ratio<=0.5)]<-0

conj_tf_change_points_old<-function(sample,tau,ratio){
  
  new.sample<-sample
  
  id<-c(which(ratio>0.5),n.grid+1)
  for (i in 1:(length(id)-1)){
    new.sample$x[new.sample$x>=grid[id[i]] & new.sample$x<(grid[id[i+1]])]<-grid[id[i]]
  }
  new.sample$x[new.sample$x==grid[id[i+1]]]
  
  sample<-new.sample
  
  grid<-unique(sample$x)
  n.grid<-length(grid)-1
  
  n01<-table(sample$x)
  idn<-which(as.character(grid) %in% names(n01))#this finds which elements are not observed
  n1<-rep(0,n.grid+1)
  n1[idn]<-n01
  yy<-tapply(sample$y,sample$x, FUN=sum)
  group.sums.y<-rep(0,n.grid+1)
  group.sums.y[idn]<-yy
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau^2*n1[n.grid+1]^2/(n1[n.grid+1]*tau^2+1) #this is going down
  for (j1 in (n.grid):1){
    new.lower<-tau^2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))^2/(tau^2*(sum(n1[j1:(n.grid+1)])-(sum(lower)))+1)
    lower<-c(new.lower,lower)
  }
  sum.lower<-rev(cumsum(rev(lower)))
  marg.weigths<-matrix(rep(c(sum.lower[3:length(sum.lower)],0),n.grid+1),ncol=n.grid+1)
  
  
  #Find the weighted means
  s<-data.frame(sample)
  local<-aggregate(s[,2],list(s$x),mean)$x
  local.means<-rep(0,n.grid+1) #CHECK IF IT IS OK TO PUT 0
  local.means[idn]<-local
  
  #Compute group local means
  lower.local.means<-local.means[n.grid+1] #this is going down
  for (j1 in (n.grid):1){
    new.lower.mean<-(sum(local.means[j1:(n.grid+1)]*n1[j1:(n.grid+1)])-sum(lower.local.means*lower[(j1+1):(n.grid+1)]))/(sum(n1[j1:(n.grid+1)])-(sum(lower[(j1+1):(n.grid+1)])))
    lower.local.means<-c(new.lower.mean,lower.local.means)
  }
  
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-rev(cumsum(rev(mean.disc)))
  marg.means<-matrix(rep(c(sum.mean.disc[3:length(sum.mean.disc)],0),n.grid+1),ncol=n.grid+1)
  
  #Define M and GAMMA and Y matrices
  sum.inv<-rev(cumsum(rev(n1)))
  sum.matrix<-matrix(rep(sum.inv,n.grid),ncol=length(sum.inv),byrow=TRUE)
  GAM<-matrix(1,nrow=n.grid,ncol=n.grid+1)
  M<-matrix(0,nrow=n.grid,ncol=n.grid+1)
  sum.y<-rev(cumsum(rev(group.sums.y)))
  Y<-matrix(rep(sum.y,n.grid),ncol=length(sum.y),byrow=TRUE)
  
  #corrected n sums
  sum.matrix.prime<-sum.matrix-marg.weigths
  Y.prime<-Y-marg.means
  
  M[,1]<-tau^2/(tau^2*sum.matrix.prime[,1]*GAM[,1]+1)
  GAM[,2]<-GAM[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]^2
  Y.prime[,2]<-Y.prime[,2]-sum.matrix.prime[,2]*M[,1]*GAM[,1]*Y.prime[,1]
  for (i in 3:(n.grid+1)){
    M[,(i-1)]<-tau^2/(tau^2*sum.matrix.prime[,(i-1)]*GAM[,(i-1)]+1)
    GAM[,i]<-GAM[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]^2)
    Y.prime[,i]<-Y.prime[,i]-sum.matrix.prime[,i]*rowSums(M[,1:(i-1)]*GAM[,1:(i-1)]*Y.prime[,1:(i-1)])
  }
  
  
  #Find posterior means
  id<-seq(2,n.grid+1)
  post.means<-(sum(n1*local.means)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))/(lower[1]+tau2^(-1))
  for (i in id){
    post.means<-c(post.means,Y.prime[i-1,i]/(sum.matrix.prime[i-1,i]*GAM[i-1,i-1]+tau2^(-1)))
  }
  
  
  #Fix the error
  post.means[which(ratio<=0.5)]<-0
  return(post.means)
  
  

  
}



solocp_single<-function(y,sigma,q=0.1,tau2=NULL,tau2.spike=NULL,tau2.slab=NULL){
  
  sigma2 <- sigma^2
  n<-length(y)
  grid <- seq(1,n)
  n.grid <- length(grid)-1
  n1<-rep(1,n.grid+1)
  
  if (is.null(tau2)){ tau2=2/sqrt(n)}
  if (is.null(tau2.spike)){tau2.spike=1/n}
  if (is.null(tau2.slab)){tau2.slab=n}
  
  weight<-matrix(1,nrow=n.grid+1,ncol=n.grid+1)
  
  #Fill the weight matrix
  lower<-tau2/(tau2+sigma2) #this is going down
  lower.local.means<-y[n.grid+1] #this is going down
  parsumLOW <- lower 
  for (j1 in (n-1):1){
    new.lower<-tau2*((n-j1+1)-parsumLOW)^2/(tau2*((n-j1+1)-(parsumLOW))+sigma2)
    lower<-c(new.lower,lower)
    parsumLOW <- parsumLOW + new.lower
    # if (abs(new.lower-1)<.000001){
    #   lower <- c(rep(1,n-length(lower)),lower)
    #   break}
  }
  revlower <- revcumsum(lower)
  revy <- revcumsum(y)
  lower.local.means<-y[n.grid+1] #this is going down
  parsum <- lower.local.means*lower[n]
  for (j1 in (n-1):1){
    new.lower.mean <- (revy[j1]-parsum)/((n-j1+1)-revlower[j1+1])
    lower.local.means<-c(new.lower.mean,lower.local.means)
    parsum <- parsum + lower[j1]*new.lower.mean
  }
  
  sum.lower<-revcumsum(lower)
  mean.disc<-lower.local.means*lower
  sum.mean.disc<-revcumsum(mean.disc)
  marg.mean <- c(sum.mean.disc[2:length(sum.mean.disc)],0)
  
  
  marg.weight <- c(sum.lower[2:length(sum.lower)],0)
  #Define M and GAMMA and Y matrices
  sum.inv<-seq(n,1)
  
  GAM<-matrix(NA,nrow=n,ncol=n)
  M<-matrix(NA,nrow=n,ncol=n)
  Y<-matrix(NA,nrow=n,ncol=n)
  
  sum.y<-rev(cumsum(rev(y)))
  
  #param
  GAM[,1] <- 1
  Y[,1] <- sum.y[1]-marg.mean
  M[,1]<-1/((sum.inv[1] - marg.weight)*GAM[,1]+sigma2*tau2^(-1))
  GAM[2:n,2]<-1-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]
  Y[2:n,2]<-sum.y[2]-marg.mean[2:n]-(sum.inv[2] - marg.weight[2:n])*M[2:n,1]*GAM[2:n,1]*Y[2:n,1]
  #weight
  w.spike<-sqrt(tau2.spike^(-1)/(lower[1]+tau2.spike^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.spike^(-1)))#maybe wrong but it does not matter
  w.slab<-sqrt(tau2.slab^(-1)/(lower[1]+tau2.slab^(-1)))*exp(1/(2*sigma2)*(sum(y)-sum(lower.local.means[2:(n.grid+1)]*lower[2:(n.grid+1)]))^2/(lower[1]+tau2.slab^(-1)))#maybe wrong but it does not matter
  parsumMGAM2 <- M[1:n,1]*GAM[1:n,1]^2
  parsumMGAMY <- M[1:n,1]*GAM[1:n,1]*Y[1:n,1]
  for (i in 3:(n-1)){
    #param
    M[(i-1):n,(i-1)]<-1/((sum.inv[(i-1)] - marg.weight[(i-1):n])*GAM[(i-1):n,(i-1)]+sigma2*tau2^(-1))
    parsumMGAM2 <- parsumMGAM2 +  M[1:n,(i-1)]*GAM[1:n,(i-1)]^2
    GAM[i:n,i]<-1-(sum.inv[i] - marg.weight[i:n])*parsumMGAM2[i:n]
    parsumMGAMY <- parsumMGAMY + M[1:n,(i-1)]*GAM[1:n,(i-1)]*Y[1:n,(i-1)]
    Y[i:n,i]<-sum.y[i]-marg.mean[i:n]-(sum.inv[i] - marg.weight[i:n])*parsumMGAMY[i:n]
    #weight
  }
  
  #param
  M[(n-1):n,(n-1)]<-1/((sum.inv[(n-1)] - marg.weight[(n-1):n])*GAM[(n-1):n,(n-1)]+sigma2*tau2^(-1))
  GAM[n,n]<-1-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]^2)
  Y[n,n]<-sum.y[i]-marg.mean[n]-(sum.inv[n] - marg.weight[n])*sum(M[n,1:(n-1)]*GAM[n,1:(n-1)]*Y[n,1:(n-1)])
  
  idGAM <- matrix(c(seq(2,n),seq(1,n-1)),ncol=2)
  w.spike<-c(w.spike,sqrt(tau2.spike^(-1)/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.spike^(-1)))*exp(1/(2*sigma2)*diag(Y)[-1]^2/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.spike^(-1))))
  w.slab<-c(w.slab,sqrt(tau2.slab^(-1)/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.slab^(-1)))*exp(1/(2*sigma2)*diag(Y)[-1]^2/((sum.inv[2:n] - marg.weight[2:n])*GAM[idGAM]+tau2.slab^(-1))))
  
  
  
  ratio<-q*w.slab/(q*w.slab+(1-q)*w.spike)
  id.inf <- which(w.slab==Inf)
  ratio[id.inf] <- 1
  
  
  return(list(w.slab=w.slab,
              w.spike=w.spike,
              ratio=ratio))
  
  
}


revcumsum <- function(x){
  return(rev(cumsum(rev(x))))
}
