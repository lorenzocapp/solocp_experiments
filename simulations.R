#### simulation trial trend filtering#####
#setwd("/Users/lorenzocappello/My Drive/Statistics/change_points/solocp_piecewise_constant/tf5_review/")
#setwd("/Users/cappello/Google Drive/postdoc Palacios/trend_filtering/code2")

source("functions_tf.R")
source("functions_server.R")
source("scenarios.R")
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
filename<-"output/"

#Parallel variable.
arg<-commandArgs(TRUE)
idx<-as.integer(arg[1])
seedz<-seq(1,100)
modelz <- seq(1,10)
comb<-expand.grid(seedz,modelz)
i<-comb[idx,1]
model <- comb[idx,2]

if (model<=4){
  jlist <- model
} else if (model==5){
  jlist <- c(5,6)
} else if (model==6){
  jlist <- c(7,8)
} else if (model > 6){
  jlist <- model +2 
}



##########################
### Start Simul Study   ##
##########################

for (j in jlist){
  print(j)
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
} else if (j==9){data <- fms.out()
data.name <-"fms.out"
} else if (j==10){data <- fms()
data.name <-"fms"
} else if (j==11){data <- fms.lap()
data.name <-"fms.lap"
} else if (j==12){data <- fms.studT()
data.name <-"fms.studT"
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


##########################
### Solo Spike and slab ##
##########################

#Parameter
tau2=2/sqrt(n)
tau2.spike=1/n
tau2.slab=n
sigma2=data$sigma^2
sigma2=sig2hat


#sigma2=mad(sample$y)^2
##########################
####         1        ####  
for (q in c(0.05,0.1,0.2,0.3)){
  for (del in c(1,2,3,4,5)){
method <- paste("solo-q",q,"-del-",del,sep="")
start_time <- Sys.time()
tf <- solocp_single(sample$y,sig2hat^(1/2),q,tau2,tau2.spike,tau2.slab)
cp <- subset_changepoints(tf$ratio,del)
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

  }
}

##########################
##########  wbs ##########
##########################
method <- "wbs_bic"
start_time <- Sys.time()
w <- wbs(sample$y)
w.cpt <- changepoints(w,penalty="bic.penalty")
cp = sort(w.cpt$cpt.ic$bic.penalty)+1
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



method <- "wbs_sic"
start_time <- Sys.time()
w <- wbs(sample$y)
w.cpt <- changepoints(w)
cp = sort(w.cpt$cpt.ic$ssic.penalty)+1
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

method <- "wbs_mbic"
start_time <- Sys.time()
w <- wbs(sample$y)
w.cpt <- changepoints(w)
cp = sort(w.cpt$cpt.ic$mbic.penalty)+1
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
cp<-which(abs(diff(fitted(smuceR(sample$y, 1:n, family="gauss",sqrt(sig2hat)))))>0)+1 #gives the change points of smuce
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
method <- "ebpiece"
start_time <- Sys.time()
o <- ebpiece(sample$y, sig2=sig2hat, 0.99, v=2*sig2hat, lambda=1, M=10000)
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


method <- "ebpiece_new"

start_time <- Sys.time()
a=fusedlasso1d(sample$y)
cv<-cv.trendfilter(a)
id1se<-which(abs(as.numeric(colnames(a$beta))-cv$lambda.1se)==min(abs(as.numeric(colnames(a$beta))-cv$lambda.1se)))
B <- rep(0,length(sample$y)-1)
B[which(diff(a$beta[,id1se])>10^-5)] <-1
o <- ebpiece_new(sample$y, sig2=sig2hat, 0.99, v=2*sig2hat, lambda=5, M=5000,B)
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
#########  pelt #########
##########################

method <- "pelt_mbic_del2"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "MBIC",minseglen = 2)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



method <- "pelt_sic_del2"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "SIC",minseglen = 2)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "pelt_bic_del2"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "BIC",minseglen = 2)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)




method <- "pelt_mbic_del3"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "BIC",minseglen = 3)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



method <- "pelt_sic_del3"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "SIC",minseglen = 3)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



method <- "pelt_bic_del3"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "BIC",minseglen = 3)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)

##########
#PELT larger spacing 
#########


method <- "pelt_mbic_del5"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "BIC",minseglen = 5)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)



method <- "pelt_sic_del5"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "SIC",minseglen = 5)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)


method <- "pelt_bic_del5"
start_time <- Sys.time()
cp <- cpt.mean(sample$y/sqrt(sig2hat), method="PELT",penalty = "BIC",minseglen = 5)@cpts+1
cp <- cp[-length(cp)]
end_time <- Sys.time()


#Evaluation
haus.dir <- dist_change_points(cp,data$change.points)
haus.inv <- dist_change_points(data$change.points,cp)
err.cp <- length(data$change.points)-length(cp)
N.vec <- N.vec_gen(est.cp=cp,true.cp=data$change.points)
N.vec.inv <- N.vec_gen.inv(est.cp=cp,true.cp=data$change.points)
sum <- c(method,data.name,seed=i,as.numeric(end_time-start_time,units="secs"),haus.dir,haus.inv,err.cp,N.vec,N.vec.inv)
df <- rbind(df,sum)




if (j>4 & j<=8){ #Skipping BASAD for large sample sizes
##########################
#########  basad #########
##########################
X <- matrix(1,ncol=n,nrow=n)
X[upper.tri(X)] <- 0



##########################
####         1        ####  


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


}

##########################
######     R-FPOP   ######
##########################


method <- "r-fpop"
start_time <- Sys.time()
res.out <- fpop_intern(sample$y/sqrt(sig2hat),  test.stat="Outlier", pen.value=2*log(length(sample$y)), 
                      lthreshold=3)
cp <- res.out$cpts[-length(res.out$cpts)]+1
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

}
