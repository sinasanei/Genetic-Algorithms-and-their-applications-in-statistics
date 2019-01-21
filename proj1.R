#install.packages("smoothmest")
#install.packages("extraDistr")
#install.packages("truncdist")
library(MASS)
#library(smoothmest)
library(stats)
library(extraDistr)

set.seed(106)
eps=runif(15)
dexp_err = -log(log(eps^-1))
f = function( x, alpha , beta, gamma){
  alpha*exp(-exp(beta-gamma*x))
}
x = seq(from=1,to=15)
y = f(x, 22.37 , 2.14 , 0.395) + dexp_err
plot(x,y,type = "l")
bs = cbind(x,y)
for (i in 1:1000){
  x_bs = sample(x,15,replace = TRUE)
  y_bs= f(x_bs, 22.37 , 2.14 , 0.395) + rgumbel(15, mu = 0, sigma = 1)
  bs = rbind(bs , cbind(x_bs,y_bs))
}
bs=data.frame(bs)
avr = NULL
low=NULL
upp=NULL
for(i in 1:15){
  bs_s = bs[which(bs$x==i),]
  a=mean(bs_s$y)
  avr[i] = a
  low[i] = a - 1.96*sd(bs_s$y)
  upp[i] = a + 1.96*sd(bs_s$y)
}
plot(x,avr,type = "l", ylim =range(-2,25), ylab = "Y")
lines(x,low, type="b", pch=22, col="red", lty=2)
lines(x,upp, type="b", pch=22, col="red", lty=2)
points(x,y,col = "blue", pch = 20)


#############
### MLE for parameter 
############
mle= function(x,y,alpha,beta,gamma){ # parameters constrained to be non-negative
  ll=0
  for ( i in 1:15){
    l= -abs(y[i]-alpha*exp(-exp(beta-gamma*x[i])))
    ll=ll+l
  }
  return(ll)
}
################
###initial population
population = cbind(a=runif(200,min=15,max=30),b=runif(200,min=0,max=5),c=runif(200,min=0,max=1) )
###########
##selection##
###########
selec = function(population ){
new_pop = matrix(NA,200,3) # initialize next generation
for( i in 1:100){
  f1= mle(x,y, population[i,1], population[i,2], population[i,3])
  f2= mle(x,y, population[i+100,1], population[i+100,2], population[i+100,3])
  if (f1 > f2){new_pop[i,]= population[i,]}
  else {new_pop[i,]= population[i+100,]}
}
return(new_pop)
}
#############
##recombination 
#############
library(truncdist)
recomb = function() { 
for(j in 1:100){
  s=sample(1:100,2)
  alphas = new_pop[s,1]
  betas =  new_pop[s,2]
  gammas =  new_pop[s,3]
  a_d = alphas[1]
  a_m = alphas[2]
  b_d = betas[1]
  b_m = betas[2]
  g_d = gammas[1]
  g_m = gammas[2]
  if(a_d > a_m){a_baby = (a_d+a_m)/2 + rtrunc(1, spec="cauchy", 
               a=(30-a_d-a_m)/(a_d-a_m) ,b=(50-a_d-a_m)/(a_d-a_m))*(a_d-a_m)/2}
  else if(a_d < a_m){a_baby = (a_d+a_m)/2 + rtrunc(1, spec="cauchy", 
               a=(50-a_d-a_m)/(a_d-a_m) ,b=(30-a_d-a_m)/(a_d-a_m))*(a_d-a_m)/2}
  else {a_baby = 0.5*(a_d+a_m) }
  if(b_d>b_m){b_baby = (b_d+b_m)/2 + rtrunc(1, spec="cauchy",
               a=(-b_d-b_m)/(b_d-b_m) ,b=(10-b_d-b_m)/(b_d-b_m))*(b_d-b_m)/2} 
  else if(b_d < b_m){b_baby = (b_d+b_m)/2 + rtrunc(1, spec="cauchy",
               a=(10-b_d-b_m)/(b_d-b_m),b=(-b_d-b_m)/(b_d-b_m))*(b_d-b_m)/2}
  else {b_baby = 0.5*(b_d+b_m) }
  if(g_d>g_m){g_baby = (g_d+g_m)/2 + rtrunc(1, spec="cauchy",
               a=(-g_d-g_m)/(g_d-g_m) ,b=(2-g_d-g_m)/(g_d-g_m))*(g_d-g_m)/2} 
  else if(g_d < g_m){g_baby = (g_d+g_m)/2 + rtrunc(1, spec="cauchy",
              a=(2-g_d-g_m)/(g_d-g_m) ,b=(-g_d-g_m)/(g_d-g_m))*(g_d-g_m)/2} 
  else {g_baby = 0.5*(g_d+g_m) }
  
  new_pop[j+100,1] = a_baby 
  new_pop[j+100,2] = b_baby
  new_pop[j+100,3] = g_baby
}
return(new_pop)
}
################
### random select 10 member from each generation and keep
################
num_cycle = 100 # 
##dim = c(10,3,num_cycle)
gen_num = 1
big_gens = list(list())

id=sample(1:200,10,replace=FALSE)
big_gens[[1]][[1]]= population[id,]


#gen_sub(population)

##############
###mutation####
###############

#############3
### Evolution
##########

for (m in 2:num_cycle){
  new_pop = selec(population)
  new_pop = recomb()
  ##########
  
  population = new_pop
  gen_num = m
  id=sample(1:200,10,replace=FALSE)
  gens[[m]]= population[id,]
  #new_pop = matrix(NA,200,3)
}
sin=NULL
for (i in 1:100){
  sin[i]= gens[[i]][5,2]
}
plot(sin,type = "l")
#####################
####################
seed = sample(1:999999999, 1000, replace = FALSE)
big_gens=list()
for (k in 1:1000){
  set.seed(seed[k])
  population = cbind(a=runif(200,min=15,max=30),b=runif(200,min=0,max=5),c=runif(200,min=0,max=1) )
  id=sample(1:200,10,replace=FALSE)
  big_gens[[k]]= list()
  big_gens[[k]][[1]]=population[id,]
  for (m in 2:num_cycle){
  
    new_pop = selec(population)
    new_pop = recomb()
    ##########
  
    population = new_pop
    gen_num = m
    id=sample(1:200,10,replace=FALSE)
    big_gens[[k]][[m]]= population[id,]
    #new_pop = matrix(NA,200,3)
  }
}
#save(big_gens, file="proj1_result.RData")
alpha=rep(NA,1000)
beta=rep(NA,1000)
gamma=rep(NA,1000)
for (i in 1:1000){
  alpha[i]=big_gens[[i]][[100]][10,1]
  beta[i]= big_gens[[i]][[100]][10,2] 
  gamma[i]=  big_gens[[i]][[100]][10,3]
}
hist(alpha,breaks = "scott",col="yellow", xlab = bquote(alpha),main="")
hist(beta,breaks = "scott",xlab = bquote(beta),main="")
hist(gamma,breaks = "scott",xlab = bquote(gamma),main="")
mean(alpha)
mean(beta)
mean(gamma)
par(mfrow=c(1,3))
hist(alpha,breaks = "scott",col="yellow",xlab = bquote(alpha),main="")
hist(beta,breaks = "scott",col="yellow",xlab = bquote(beta),main="")
hist(gamma,breaks = "scott",col="yellow",xlab = bquote(gamma),main="")
