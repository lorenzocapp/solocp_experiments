#### simulation trial trend filtering#####

source("functions_server_section_4_1.R")
source("functions_tf_section_4_1.R")
source("simulation_scenarios_section_4_1.R")
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
model <- seq(9,12)
comb<-expand.grid(seedz,model)
i<-comb[idx,1]
j<-comb[idx,2]




##########################
### Start Simul Study   ##
##########################


#Load Data
set.seed(i)
if (j==1) {data <- blocks7.out()
data.name <-"blocks7.out"
} else if (j==2){data <- blocks7()
data.name <-"blocks7"
} else if (j==3){data <- blocks7.lap()
data.name <-"blocks7.lap"
} else if (j==4){data <- blocks7.studT()
data.name <-"blocks7.studT"
} else if (j==5){data <- smallwideteeth.25.out()
data.name <-"smallwideteeth.25.out"
} else if (j==6){data <- smallwideteeth.25()
data.name <-"smallwideteeth.25"
} else if (j==7){data <- smallwideteeth.25.lap()
data.name <-"smallwideteeth.25.lap"
} else if (j==8){data <- smallwideteeth.25.studT()
data.name <-"smallwideteeth.25.studT"
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
df <- rep(NA,21)

sig2hat
##########################
### Solo Spike and slab ##
##########################

#Parameter
tau2=2/sqrt(n)
tau2.spike=1/n
tau2.slab=n
#sigma2=data$sigma^2
sigma2=sig2hat


##########################
####         1        ####  
q=0.05
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=5)
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
####         2        ####  
q=0.1
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=5)
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
####         3        ####  
q=0.2
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=5)
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
####         4        ####  
q=0.3
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=5)
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
####         5        ####  
q=0.5
method <- paste("solo-q",q,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=5)
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
##########  wbs ##########
##########################
method <- "wbs"
start_time <- Sys.time()
w <- wbs(sample$y)
w.cpt <- changepoints(w,penalty="bic.penalty")
cp = sort( w.cpt$cpt.ic$bic.penalty)+1
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
#########  smuce #########
##########################
method <- "smuce"
start_time <- Sys.time()
cp<-which(abs(diff(fitted(smuceR(sample$y, 1:n, family="gauss"))))>0)+1 #gives the change points of smuce
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


############################
#########  ebpiece #########
############################
method <- "ebpiece_new_2"

start_time <- Sys.time()
B <- rep(0,length(sample$y)-1)
B[which(diff(a$beta[,id1se])>10^-5)] <-1
o <- ebpiece_new(sample$y, sig2=sig2hat, 0.99, v=2*sig2hat, lambda=2, M=10000,B)
cp <- which(diff(apply(o$B, 2, mean))>0)+2
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "ebpiece_new_5"

start_time <- Sys.time()
B <- rep(0,length(sample$y)-1)
B[which(diff(a$beta[,id1se])>10^-5)] <-1
o <- ebpiece_new(sample$y, sig2=sig2hat, 0.99, v=2*sig2hat, lambda=5, M=10000,B)
cp <- which(diff(apply(o$B, 2, mean))>0)+2
end_time <- Sys.time()

#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

method <- "ebpiece_new_10"

start_time <- Sys.time()
B <- rep(0,length(sample$y)-1)
B[which(diff(a$beta[,id1se])>10^-5)] <-1
o <- ebpiece_new(sample$y, sig2=sig2hat, 0.99, v=2*sig2hat, lambda=10, M=10000,B)
cp <- which(diff(apply(o$B, 2, mean))>0)+2
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
######     R-FPOP   ######
##########################


method <- "r-fpop"
start_time <- Sys.time()
res.l2 <- Rob_seg.std(x = sample$y/sqrt(sig2hat), loss = "Outlier",  lambda=2*log(length(sample$y)),lthreshold=3*sqrt(sig2hat))
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

est.sd <- mad(diff(sample$y)/sqrt(2))
method <- "r-fpop_mad"
start_time <- Sys.time()
res.l2 <- Rob_seg.std(x = sample$y/est.sd, loss = "Outlier",  lambda=2*log(length(sample$y))lthreshold=3*est.sd)
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
#########  pelt #########
##########################

method <- "pelt"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/mad(diff(sample$y)/sqrt(2)), method="PELT")@cpts
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
X <- matrix(1,ncol=n,nrow=n)
X[upper.tri(X)] <- 0



##########################
####         1        ####  
K = floor(.05*n)
method <- paste("basad-K0.05",sep="")

start_time <- Sys.time()
Z0 = array(0,(n))
bsd = withTimeout(basad(X,sample$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
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
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         2        ####  
K = floor(.1*n)
method <- paste("basad-K0.1",sep="")

start_time <- Sys.time()
Z0 = array(0,(n))
bsd = withTimeout(basad(X,sample$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
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
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         3        ####  
K = floor(.2*n)
method <- paste("basad-K0.2",sep="")

start_time <- Sys.time()
Z0 = array(0,(n))
bsd = withTimeout(basad(X,sample$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
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
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


##########################
####         4        ####  
K = floor(.3*n)
method <- paste("basad-K0.3",sep="")

start_time <- Sys.time()
Z0 = array(0,(n))
bsd = withTimeout(basad(X,sample$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
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
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########################
####         5        ####  
K = floor(.5*n)
method <- paste("basad-K0.5",sep="")

start_time <- Sys.time()
Z0 = array(0,(n))
bsd = withTimeout(basad(X,sample$y,Z0,returnB = 1,K = K,del = 0.1,sig = sig2hat),timeout=7200,onTimeout = "silent")
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
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



####################################
### Solo Spike and slab TEST DEL  ##
###################################

#Parameter
tau2=2/sqrt(n)
tau2.spike=1/n
tau2.slab=n
#sigma2=data$sigma^2
sigma2=sig2hat


##########################
####         1        ####  
q=0.1
del =1
method <- paste("solo-q",q,"-del",del,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=del)
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
####         2        ####  
q=0.1
del =3
method <- paste("solo-q",q,"-del",del,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=del)
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
####         3        ####  
q=0.1
del =7
method <- paste("solo-q",q,"-del",del,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=del)
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
####         4        ####  
q=0.1
del =9
method <- paste("solo-q",q,"-del",del,sep="")
start_time <- Sys.time()
tf<-solo_spike_slab_grid_free(sample,tau2,tau2.spike,tau2.slab,q,sigma2)
cp <- subset_changepoints(tf$ratio,del=del)
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
######     R-FPOP   ######
##########################


method <- "r-fpop"
start_time <- Sys.time()
res.l2 <- Rob_seg.std(x = sample$y/sqrt(sig2hat), loss = "L2",  lambda=2*log(length(sample$y)))
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



#################
### SAVE ########
#################
df <- df[-1,]
outputName=paste(filename,"out-seed",i,"-j-",j,"-v1.RData",sep="")
outputPath=file.path(outputName)
save(df,file=outputPath)


