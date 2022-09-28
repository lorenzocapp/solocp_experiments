#### simulation trial trend filtering#####

source("functions_tf_binned.R")
source("functions_server_binned.R")
source("scenarios_binned.R")
library(genlasso)
library(wbs)
library(stepR)
library(changepoint)
library(cumSeg)
library(strucchange)
library(R.utils)
library(rmutil)
library(robseg)




#File to save (create a corresponding folder in scratch)
filename<-""

#Parallel variable.
arg<-commandArgs(TRUE)
idx<-as.integer(arg[1])
seedz<-seq(1,100)
model <- seq(1,4)
comb<-expand.grid(seedz,model)
#idx<-1
i<-comb[idx,1]
j<-comb[idx,2]




##########################
### Start Simul Study   ##
##########################


#Load Data
set.seed(i)
if (j==1) {data <- blocks7.out.mod.r()
data.name <-"blocks7.out"
} else if (j==2){data <- blocks7.r()
data.name <-"blocks7"
} else if (j==3){data <- blocks7.lap.r()
data.name <-"blocks7.lap"
} else if (j==4){data <- blocks7.studT.r()
data.name <-"blocks7.studT"
}

sample <-data$sample
n<-dim(data$sample)[1]

#Compute sigma hat squared
a=fusedlasso1d(sample$y)
cv<-cv.trendfilter(a)
id<-which(abs(as.numeric(colnames(a$beta))-cv$lambda.min)==min(abs(as.numeric(colnames(a$beta))-cv$lambda.min)))
active <- sum(diff(a$beta[,id])>10^-3)
sig2hat <- sum((sample$y-a$beta[,id])^2)/(n-active)
#Create data frame 
df <- rep(NA,14)


n.grid <- 200
data$change.points.grid <- cp_grid(data$change.points,dim(sample)[1],n.grid)$cp.grid

##########################
### Solo Spike and slab ##
##########################

#Parameter
tau2=2/sqrt(n.grid)
tau2.spike=1/n.grid
tau2.slab=n.grid
#sigma2=data$sigma^2
sigma2=sig2hat


##########################
####         1        ####  
q=0.05
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_brack(sample,tau2,tau2.spike,tau2.slab,q,sigma2,n.grid)
cp <- subset_changepoints(tf$ratio,del=5)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         2        ####  
q=0.1
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_brack(sample,tau2,tau2.spike,tau2.slab,q,sigma2,n.grid)
cp <- subset_changepoints(tf$ratio,del=5)
cp
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         3        ####  
q=0.2
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_brack(sample,tau2,tau2.spike,tau2.slab,q,sigma2,n.grid)
cp <- subset_changepoints(tf$ratio,del=5)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         4        ####  
q=0.3
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_brack(sample,tau2,tau2.spike,tau2.slab,q,sigma2,n.grid)
cp <- subset_changepoints(tf$ratio,del=5)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         5        ####  
q=0.5
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_brack(sample,tau2,tau2.spike,tau2.slab,q,sigma2,n.grid)
cp <- subset_changepoints(tf$ratio,del=5)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
##########  wbs ##########
##########################
method <- "wbs_loc.means"
start_time <- Sys.time()
loc.means <- local_means_brack(sample,n.grid)
w <- wbs(loc.means)
w.cpt <- changepoints(w,penalty="bic.penalty")
cp = sort( w.cpt$cpt.ic$bic.penalty)+1
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "wbs_post"
start_time <- Sys.time()
sample2 <- sample[order(sample$x),]
w <- wbs(sample2$y)
w.cpt <- changepoints(w,penalty="bic.penalty")
cp0 = sort( w.cpt$cpt.ic$bic.penalty)+1
cp <- unique(cp_grid(cp0,dim(sample)[1],n.grid)$cp.grid)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
#########  smuce #########
##########################
method <- "smuce_loc.means"
start_time <- Sys.time()
cp<-which(abs(diff(fitted(smuceR(loc.means, 1:n.grid, family="gauss"))))>0)+1 #gives the change points of smuce
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "smuce_post"
start_time <- Sys.time()
cp0<-which(abs(diff(fitted(smuceR(sample2$y, 1:n, family="gauss"))))>0)+1 #gives the change points of smuce
cp <- unique(cp_grid(cp0,dim(sample)[1],n.grid)$cp.grid)
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

############################
#########  ebpiece #########
############################
method <- "ebpiece_loc.means"
start_time <- Sys.time()
o <- ebpiece(loc.means, sig2=sig2hat, 0.99, v=1, lambda=2, M=10000)
cp <- which(diff(apply(o$B, 2, mean))>0)+2
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "ebpiece_post"
start_time <- Sys.time()
o <- ebpiece(sample2$y, sig2=sig2hat, 0.99, v=1, lambda=2, M=10000)
cp0 <- which(diff(apply(o$B, 2, mean))>0)+2
cp <- unique(cp_grid(cp0,dim(sample)[1],n.grid)$cp.grid)
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
#########  pelt #########
##########################

method <- "pelt_loc.means"
start_time <- Sys.time()
cp <- cpt.mean(loc.means/mad(diff(loc.means)/sqrt(2)), method="PELT")@cpts
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "pelt_post"
start_time <- Sys.time()
cp0 <- cpt.mean(sample2$y/mad(diff(sample2$y)/sqrt(2)), method="PELT")@cpts
cp <- unique(cp_grid(cp0,dim(sample)[1],n.grid)$cp.grid)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
######     R-FPOP   ######
##########################


method <- "r-fpop_loc.means"
start_time <- Sys.time()
res.l2 <- Rob_seg.std(x = loc.means/sqrt(sig2hat), loss = "Outlier",  lambda=2*log(length(loc.means)),lthreshold=3*sqrt(sig2hat))
cp <- res.l2$t.est[-length(res.l2$t.est)]+1
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


est.sd <- mad(diff(loc.means)/sqrt(2))
method <- "r-fpop_loc.means_mad"
start_time <- Sys.time()
res.l2 <- Rob_seg.std(x = loc.means/est.sd, loss = "Outlier",  lambda=2*log(length(loc.means)),lthreshold=3*est.sd)
cp <- res.l2$t.est[-length(res.l2$t.est)]+1
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
#########  basad #########
##########################

X <- X_mat_brack(n.grid,sample,n)



##########################
####         1        ####  
K = floor(.05*n.grid)
method <- paste("basad-K0.05",sep="")

start_time <- Sys.time()
Z0 = array(0,(n.grid))
bsd = withTimeout(basad(X,sample2$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
if (length(bsd)>0){
  Z = bsd$Z
  ZZ1 = apply(Z,2,mean)
  cp <- subset_changepoints(ZZ1,del=5)
  end_time <- Sys.time()
} else {
  cp <- NULL
  end_time <- Sys.time()
}


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         2        ####  
K = floor(.1*n.grid)
method <- paste("basad-K0.1",sep="")

start_time <- Sys.time()
Z0 = array(0,(n.grid))
bsd = withTimeout(basad(X,sample2$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
if (length(bsd)>0){
  Z = bsd$Z
  ZZ1 = apply(Z,2,mean)
  cp <- subset_changepoints(ZZ1,del=5)
  end_time <- Sys.time()
} else {
  cp <- NULL
  end_time <- Sys.time()
}


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
####         3        ####  
K = floor(.2*n)
method <- paste("basad-K0.2",sep="")

start_time <- Sys.time()
Z0 = array(0,(n.grid))
bsd = withTimeout(basad(X,sample2$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
if (length(bsd)>0){
  Z = bsd$Z
  ZZ1 = apply(Z,2,mean)
  cp <- subset_changepoints(ZZ1,del=5)
  end_time <- Sys.time()
} else {
  cp <- NULL
  end_time <- Sys.time()
}


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



##########################
####         4        ####  
K = floor(.3*n)
method <- paste("basad-K0.3",sep="")

start_time <- Sys.time()
Z0 = array(0,(n.grid))
bsd = withTimeout(basad(X,sample2$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
if (length(bsd)>0){
  Z = bsd$Z
  ZZ1 = apply(Z,2,mean)
  cp <- subset_changepoints(ZZ1,del=5)
  end_time <- Sys.time()
} else {
  cp <- NULL
  end_time <- Sys.time()
}


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
####         5        ####  
K = floor(.5*n.grid)
method <- paste("basad-K0.5",sep="")

start_time <- Sys.time()
Z0 = array(0,(n.grid))
bsd = withTimeout(basad(X,sample2$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
if (length(bsd)>0){
  Z = bsd$Z
  ZZ1 = apply(Z,2,mean)
  cp <- subset_changepoints(ZZ1,del=5)
  end_time <- Sys.time()
} else {
  cp <- NULL
  end_time <- Sys.time()
}


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points.grid)
haus.inv <- dist_change_points(data$change.points.grid,cp)
err.cp <- length(data$change.points.grid)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points.grid)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points.grid)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)




#################
### SAVE ########
#################
df <- df[-1,]
outputName=paste(filename,"out-seed",i,"-j-",j,"-v1.RData",sep="")
outputPath=file.path(outputName)
save(df,file=outputPath)


