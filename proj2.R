N = 200 # population (number of chains)
d= 100 # dimention 
#initial = runif(N*d,min = 9.9 , max= 10)
cov = diag(1:d)
cov
for (i in 1:d){
  for (j in 1:d){
    if (i != j){
      cov[i,j]=0.5
    }
  }
}
cov_inv=solve(cov)
################
###multivariate normal log-lik
##############
log_lik =function(x){ ## fitness function
 #-0.5*(log(det(cov))+(t(x)%*%cov_inv%*%x)+length(x)*log(2*pi))
  -log(det(cov))-(t(x)%*%cov_inv%*%x)
}
############
###optimal factor gamma
#########3
C= 2.38 / sqrt(2*d)
seq=1:N
#########
##initial sample
###########
X = matrix(runif(N*d,min = 9.9 , max= 10),N,d) 
###########
generations = 5000 # chain length 
samples = list()
samples[[1]]= X
X_mat = X
set.seed(12)
 for ( i in 2:generations) { 
   for ( j in 1:N){
     a=sample(seq[which(seq!=j)],2,replace = FALSE )
     R1 = X_mat[a[1],]
     R2= X_mat[a[2],]
     X_p=NULL #proposal
     X_p=X_mat[j,]+ C*(R1-R2) + runif(d,min = -0.0001, max = 0.0001)
     #r = log_lik(X_p)/ log_lik(X_mat[j,])
     r = min(log_lik(X_p)- log_lik(X_mat[j,]),0)
     #selection process: accept if metropolis ration > uniform(0,1)
     if (r> log(runif(1))) { X_mat[j,] = X_p }       #(log(r) > log(runif(1))) { X_mat[j,] = X_p }
   }
   samples[[i]]=X_mat
 }

s=rep(NA,5000)
for(p in 1:5000){
  x= samples[[p]]
  s[p]=x[10,38]
}
plot(s,type = "l")  
X1 = matrix(NA,200,5000)
X2 =  matrix(NA,200,5000)
for( p in 1:5000){
  X_sam = samples[[p]]
  X1[,p]= X_sam[,1]
  X2[,p]= X_sam[,100]
  
}
mean=NULL
sd1=NULL
sd2=NULL
cov=NULL
for (s in 1:5000){
  mean[s]=mean(X1[,s])
  sd1[s]=sd(X1[,s])
  sd2[s]=sd(X2[,s])
  cov[s]= cov(X1[,s],X2[,s]) 
    
}
data1 = cbind(mean,sd1,sd2,cov)
save(data1, file="proj2_result1.RData")
####################################RUN 2222222222222222222222222222
N = 200 # population (number of chains)
d= 100 # dimention 
#initial = runif(N*d,min = 9.9 , max= 10)
cov = diag(1:d)
cov
for (i in 1:d){
  for (j in 1:d){
    if (i != j){
      cov[i,j]=0.5
    }
  }
}
cov_inv=solve(cov)
##initial sample
###########
X = matrix(runif(N*d,min = -5 , max= 15),N,d) 
###########
generations = 5000 # chain length 
samples = list()
samples[[1]]= X
X_mat = X
set.seed(12)
for ( i in 2:generations) { 
  for ( j in 1:N){
    a=sample(seq[which(seq!=j)],2,replace = FALSE )
    R1 = X_mat[a[1],]
    R2= X_mat[a[2],]
    X_p=NULL #proposal
    X_p=X_mat[j,]+ C*(R1-R2) + runif(d,min = -0.0001, max = 0.0001)
    #r = log_lik(X_p)/ log_lik(X_mat[j,])
    r = min(log_lik(X_p)- log_lik(X_mat[j,]),0)
    #selection process: accept if metropolis ration > uniform(0,1)
    if (r> log(runif(1))) { X_mat[j,] = X_p }       #(log(r) > log(runif(1))) { X_mat[j,] = X_p }
  }
  samples[[i]]=X_mat
}
X1 = matrix(NA,200,5000)
X2 =  matrix(NA,200,5000)
for( p in 1:5000){
  X_sam = samples[[p]]
  X1[,p]= X_sam[,1]
  X2[,p]= X_sam[,100]
  
}
mean2=NULL
sd12=NULL
sd22=NULL
cov2=NULL
for (s in 1:5000){
  mean2[s]=mean(X1[,s])
  sd12[s]=sd(X1[,s])
  sd22[s]=sd(X2[,s])
  cov2[s]= cov(X1[,s],X2[,s]) 
  
}
data2 = cbind(mean2,sd12,sd22,cov2)
save(data2, file="proj2_result2.RData")
#########################################RUN33333333333333333333333333333333
N = 101 # population (number of chains)
d= 100 # dimention 
#initial = runif(N*d,min = 9.9 , max= 10)
cov = diag(1:d)
cov
for (i in 1:d){
  for (j in 1:d){
    if (i != j){
      cov[i,j]=0.5
    }
  }
}
cov_inv=solve(cov)
################
###multivariate normal log-lik
##############
log_lik =function(x){ ## fitness function
  #-0.5*(log(det(cov))+(t(x)%*%cov_inv%*%x)+length(x)*log(2*pi))
  -log(det(cov))-(t(x)%*%cov_inv%*%x)
}
############
###optimal factor gamma
#########3
C= 2.38 / sqrt(2*d)
seq=1:N
#########
##initial sample
###########
X = matrix(runif(N*d,min = -5 , max= 15),N,d) 
###########
generations = 10000 # chain length 
samples = list()
samples[[1]]= X
X_mat = X
set.seed(12)
for ( i in 2:generations) { 
  for ( j in 1:N){
    a=sample(seq[which(seq!=j)],2,replace = FALSE )
    R1 = X_mat[a[1],]
    R2= X_mat[a[2],]
    X_p=NULL #proposal
    X_p=X_mat[j,]+ C*(R1-R2) + runif(d,min = -0.0001, max = 0.0001)
    #r = log_lik(X_p)/ log_lik(X_mat[j,])
    r = min(log_lik(X_p)- log_lik(X_mat[j,]),0)
    #selection process: accept if metropolis ration > uniform(0,1)
    if (r> log(runif(1))) { X_mat[j,] = X_p }       #(log(r) > log(runif(1))) { X_mat[j,] = X_p }
  }
  samples[[i]]=X_mat
}
X1 = matrix(NA,101,10000)
X2 =  matrix(NA,101,10000)
for( p in 1:10000){
  X_sam = samples[[p]]
  X1[,p]= X_sam[,1]
  X2[,p]= X_sam[,100]
  
}
mean3=NULL
sd13=NULL
sd23=NULL
cov3=NULL
for (s in 1:10000){
  mean3[s]=mean(X1[,s])
  sd13[s]=sd(X1[,s])
  sd23[s]=sd(X2[,s])
  cov3[s]= cov(X1[,s],X2[,s]) 
  
}
data3 = cbind(mean3,sd13,sd23,cov3)
save(data3, file="proj2_result3.RData")
##########################RUN444444444444444444444444444
N = 1000 # population (number of chains)
d= 100 # dimention 
#initial = runif(N*d,min = 9.9 , max= 10)
cov = diag(1:d)
cov
for (i in 1:d){
  for (j in 1:d){
    if (i != j){
      cov[i,j]=0.5
    }
  }
}
cov_inv=solve(cov)
################
###multivariate normal log-lik
##############
log_lik =function(x){ ## fitness function
  #-0.5*(log(det(cov))+(t(x)%*%cov_inv%*%x)+length(x)*log(2*pi))
  -log(det(cov))-(t(x)%*%cov_inv%*%x)
}
############
###optimal factor gamma
#########3
C= 2.38 / sqrt(2*d)
seq=1:N
#########
##initial sample
###########
X = matrix(runif(N*d,min = -5 , max= 15),N,d) 
###########
generations = 1000 # chain length 
samples = list()
samples[[1]]= X
X_mat = X
set.seed(12)
for ( i in 2:generations) { 
  for ( j in 1:N){
    a=sample(seq[which(seq!=j)],2,replace = FALSE )
    R1 = X_mat[a[1],]
    R2= X_mat[a[2],]
    X_p=NULL #proposal
    X_p=X_mat[j,]+ C*(R1-R2) + runif(d,min = -0.0001, max = 0.0001)
    #r = log_lik(X_p)/ log_lik(X_mat[j,])
    r = min(log_lik(X_p)- log_lik(X_mat[j,]),0)
    #selection process: accept if metropolis ration > uniform(0,1)
    if (r> log(runif(1))) { X_mat[j,] = X_p }       #(log(r) > log(runif(1))) { X_mat[j,] = X_p }
  }
  samples[[i]]=X_mat
}
X1 = matrix(NA,1000,1000)
X2 =  matrix(NA,1000,1000)
for( p in 1:1000){
  X_sam = samples[[p]]
  X1[,p]= X_sam[,1]
  X2[,p]= X_sam[,100]
  
}
mean4=NULL
sd14=NULL
sd24=NULL
cov4=NULL
for (s in 1:1000){
  mean4[s]=mean(X1[,s])
  sd14[s]=sd(X1[,s])
  sd24[s]=sd(X2[,s])
  cov4[s]= cov(X1[,s],X2[,s]) 
  
}
data4 = cbind(mean4,sd14,sd24,cov4)
save(data4, file="proj2_result4.RData")
#######################
load(file="proj2_result1.RData")
load(file="proj2_result2.RData")
load(file="proj2_result3.RData")
load(file="proj2_result4.RData")
##################
par(mfrow=c(2,2))
plot(data1[,1],xlab = "Generations",ylab = "", type = "l", ylim = c(-2,12), main = "N=200 d=100, init X=[9.9,10]")
points(data1[,2], type = "l" ,col="green")
points(data1[,3], type = "l", col="red" )
points(data1[,4], type = "l",col="blue")
plot(data2[,1],xlab = "Generations",ylab = "", type = "l", ylim = c(-2,12), main = "N=200 d=100, init X=[-5,15]")
points(data2[,2], type = "l" ,col="green")
points(data2[,3], type = "l", col="red" )
points(data2[,4], type = "l",col="blue")
plot(data3[,1],xlab = "Generations",ylab = "", type = "l", ylim = c(-2,12), main = "N=101 d=100, init X=[-5,15]")
points(data3[,2], type = "l" ,col="green")
points(data3[,3], type = "l", col="red" )
points(data3[,4], type = "l",col="blue")
plot(data4[,1],xlab = "Generations",ylab = "", type = "l", ylim = c(-2,12), main = "N=1000 d=100, init X=[-5,15]")
points(data4[,2], type = "l" ,col="green")
points(data4[,3], type = "l", col="red" )
points(data4[,4], type = "l",col="blue")
##################################
##################################
library(mcmcse)
mcse(data1[,1], size = "sqroot", g = NULL,
      method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#1
mcse(data1[,2], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#2
mcse(data1[,3], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#3
mcse(data1[,4], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#4
mcse(data2[,1], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#5
mcse(data2[,2], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#6
mcse(data2[,3], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#7
mcse(data2[,4], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#8
mcse(data3[,1], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#9
mcse(data3[,2], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#10
mcse(data3[,3], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#11
mcse(data3[,4], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#12
mcse(data4[,1], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#13
mcse(data4[,2], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#14
mcse(data4[,3], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#15
mcse(data4[,4], size = "sqroot", g = NULL,
     method = c("bm", "obm", "tukey", "bartlett"),warn = FALSE)#16
#mcerror_isadj <- mcse.initseq(x = data1, g = NULL,
 #                             level = .95, adjust = TRUE)
#mcerror_bm$cov
