#install.packages("MCMCpack")
library(MCMCpack)

est_method<-c(rep("1", 6), rep("2", 6), rep("3", 6),rep("4", 6))
##Instructor<-c(rep(c("1","1","2","2"), 3))
observations<-c(0.34,0.12,1.23,0.70,1.75,0.12,
                0.91,2.94,2.14,2.36,2.86,4.55,
                6.31,8.37,9.75,6.09,9.82,7.24,
                17.15,11.82,10.95,17.20,14.35,16.82)
df<-data.frame(est_method=as.factor(est_method), observations = observations)
I = 6
J = 4
d=7
N=  2*d #10*d
#############
alpha = 1
beta= 10
a = 1.85
b= 0.1
##########
###  p(mu,log(sigma^2),log(tau^2),theta) #####
########

loglik_pi = function(vec){ ####vec=(mu,log(sigma^2),log(tau^2),theta) 
  mu= vec[1]
  l_sigma_2=vec[2]
  l_tau_2=vec[3]
  theta = vec[4:7]
  llN_theta =0
  Z = 0
  for ( j in 1:J){
    llN_theta = llN_theta - 0.5*log(2*pi*exp(l_tau_2)) - 0.5*(theta[j]-mu)^2*(exp(l_tau_2)^-1)
    llN_y =0
    for ( i in 1:I){
      llN_y= llN_y -  0.5*log(2*pi*exp(l_sigma_2)) - 0.5*(y[i,j]-theta[j])^2*(exp(l_sigma_2)^-1) 
    }
    Z = Z + llN_theta + llN_y
    N_theta = 0
  }
  return( -2*alpha*log(sqrt(exp(l_sigma_2)))-beta*(exp(l_sigma_2)^-1)
          - 2*a*log(sqrt(exp(l_tau_2)))- b*(exp(l_tau_2)^-1) + Z )
}

##################
### initial population ##### X=()
#################
X= matrix(NA,N,d,
      dimnames = list(NULL,c("mu","log(sigma^2)","log(tau^2)","theta1","theta2","theta3","theta4")))
X[,1]=runif(N,min = -20,max = 20 )         ### mu
X[,2]=log(rinvgamma(N,alpha,scale=beta))  ### log(sigma^2)
X[,3]=log(rinvgamma(N,a, scale=b))       ### log(tau^2)
y=matrix(df$observations,6,4)
for( i in 1:N){
  X[i,4]=rnorm(1,X[i,1], sd=sqrt(exp(X[i,3])))  ## theta_1
  X[i,5]=rnorm(1,X[i,1], sd=sqrt(exp(X[i,3])))  ## theta_2
  X[i,6]=rnorm(1,X[i,1], sd=sqrt(exp(X[i,3])))  ## theta_3
  X[i,7]=rnorm(1,X[i,1], sd=sqrt(exp(X[i,3])))  ## theta_4
}

###############3
generations = 1000000 #10^6 # chain length 
samples = list()
samples[[1]]= X
X_mat = X
seq=1:N
set.seed(142)
for ( i in 2:generations) { 
  if ( i%%10==0){ 
    C= 1 
    } else {C= 2.38 / sqrt(2*d)}
  for ( j in 1:N){
    s=sample(seq[which(seq!=j)],2,replace = FALSE )
    P1 = X_mat[s[1],]
    P2= X_mat[s[2],]
    X_p=NULL #proposal
    X_p=X_mat[j,]+ C*(P1-P2) + runif(d,min = -0.0001, max = 0.0001)
    #X_p[4]=rnorm(1,X_p[1], sd=sqrt(exp(X_p[3]))) 
    #X_p[5]=rnorm(1,X_p[1], sd=sqrt(exp(X_p[3]))) 
    #X_p[6]=rnorm(1,X_p[1], sd=sqrt(exp(X_p[3]))) 
    #X_p[7]=rnorm(1,X_p[1], sd=sqrt(exp(X_p[3]))) 
    r = min( loglik_pi(X_p)-loglik_pi(X_mat[j,]),0) 
    #selection process: accept if metropolis ration > uniform(0,1)
    if (r > log(runif(1))) { X_mat[j,] = X_p }       #(log(r) > log(runif(1))) { X_mat[j,] = X_p }
  }
  samples[[i]]=X_mat
}

#######
save(samples, file="proj3_result.RData")
#########
zeta = NULL
zeta2=NULL
phi = NULL
e=sample(1:14,1,replace = FALSE)
for(k in 1:10000){
  zeta[k]= samples[[k]][e,2]
  zeta2[k]= samples[[k]][e,3]
  #zeta[k]=exp(samples[[k+100000]][e,2])/exp(samples[[k+100000]][e,3])
  #phi[k] = exp(samples[[k+100000]][e,2])/ (I*exp(samples[[k+100000]][e,3])+exp(samples[[k+100000]][e,2]))
}
plot(zeta,type = "l")
###################
load(file="proj3_result.RData")
mat2= matrix(NA,900000*14,2)

for ( i in 100001:1000000){
  ss=samples[[i]]
  mat2[(i-100001)*14-13,1]=log( exp(ss[1,2])/exp(ss[1,3]))
  mat2[(i-100001)*14-13,2]=exp(ss[1,2])/(6* exp(ss[1,3])+exp(ss[1,2]))
  mat2[(i-100001)*14-12,1]=log( exp(ss[2,2])/exp(ss[2,3]))
  mat2[(i-100001)*14-12,2]=exp(ss[2,2])/(6* exp(ss[2,3])+exp(ss[2,2]))
  mat2[(i-100001)*14-11,1]=log( exp(ss[3,2])/exp(ss[3,3]))
  mat2[(i-100001)*14-11,2]=exp(ss[3,2])/(6* exp(ss[3,3])+exp(ss[3,2]))
  mat2[(i-100001)*14-10,1]=log( exp(ss[4,2])/exp(ss[4,3]))
  mat2[(i-100001)*14-10,2]=exp(ss[4,2])/(6* exp(ss[4,3])+exp(ss[4,2]))
  mat2[(i-100001)*14-9,1]=log (exp(ss[5,2])/exp(ss[5,3]))
  mat2[(i-100001)*14-9,2]=exp(ss[5,2])/(6* exp(ss[5,3])+exp(ss[5,2]))
  mat2[(i-100001)*14-8,1]=log( exp(ss[6,2])/exp(ss[6,3]))
  mat2[(i-100001)*14-8,2]=exp(ss[6,2])/(6* exp(ss[6,3])+exp(ss[6,2]))
  mat2[(i-100001)*14-7,1]=log( exp(ss[7,2])/exp(ss[7,3]))
  mat2[(i-100001)*14-7,2]=exp(ss[7,2])/(6* exp(ss[7,3])+exp(ss[7,2]))
  mat2[(i-100001)*14-6,1]=log( exp(ss[8,2])/exp(ss[8,3]))
  mat2[(i-100001)*14-6,2]=exp(ss[8,2])/(6* exp(ss[8,3])+exp(ss[8,2]))
  mat2[(i-100001)*14-5,1]=log( exp(ss[9,2])/exp(ss[9,3]))
  mat2[(i-100001)*14-5,2]=exp(ss[9,2])/(6* exp(ss[9,3])+exp(ss[9,2]))
  mat2[(i-100001)*14-4,1]=log( exp(ss[10,2])/exp(ss[10,3]))
  mat2[(i-100001)*14-4,2]=exp(ss[10,2])/(6* exp(ss[10,3])+exp(ss[10,2]))
  mat2[(i-100001)*14-3,1]=log( exp(ss[11,2])/exp(ss[11,3]))
  mat2[(i-100001)*14-3,2]=exp(ss[11,2])/(6* exp(ss[11,3])+exp(ss[11,2]))
  mat2[(i-100001)*14-2,1]=log( exp(ss[12,2])/exp(ss[12,3]))
  mat2[(i-100001)*14-2,2]=exp(ss[12,2])/(6* exp(ss[12,3])+exp(ss[12,2]))
  mat2[(i-100001)*14-1,1]=log( exp(ss[13,2])/exp(ss[13,3]))
  mat2[(i-100001)*14-1,2]=exp(ss[13,2])/(6* exp(ss[13,3])+exp(ss[13,2]))
  mat2[(i-100001)*14,1]=log( exp(ss[14,2])/exp(ss[14,3]))
  mat2[(i-100001)*14,2]=exp(ss[14,2])/(6* exp(ss[14,3])+exp(ss[14,2]))
}
