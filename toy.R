library(binaryLogic)
library(bit)
set.seed(5)
###########################################
N= 5  #number of initial population size  #
M = 20 #10 # #number of generations  ##########
###########################################
#randomly generated initial population of size N 
################################################
populations = matrix(NA,M,N)
fitness= matrix(NA,M,N)
##############################
init=sample(2^6-1,N)
pop_bin=list()
for (p in 1:N){
  pop_bin[[p]] = as.numeric(unclass(as.binary(init[p],n=6)))
}
#pop_bin = as.binary(init,n=6)
populations[1,] = init
############################
#objective function value##
###########################
f= function(x){
  val = x^2 -42*x +152
  return(val)
}
##########################################################
## Fitness value for population members saved in a MAtrix##
##########################################################
#fit = function(population, func){
#   fitness=rep(NA,length(population))
#  for(i in 1:length(population)){
#    fitness[i] = -func(as.numeric(population[[i]]))
#    }
#  return(fitness)
#}
fitness[1,]=-f(populations[1,])

##################
######binary to numeric
##################
bintonum = function(alist){
  l = length(alist)
  row_in_pop =rep(NA,length(alist))
  for (k in 1:l){
    obj=as.numeric(alist[[k]])
    row_in_pop[k]= obj[1]*(2^5)+obj[2]*(2^4)+obj[3]*(2^3)+obj[4]*(2^2)+obj[5]*2+obj[6]
  }
  return(row_in_pop)
}
########################################
#next.gen = list() #retain fittest chromosome 
#next.gen[1]= pop[which.max(fitness[1,])]
###selection Tournament
selection = function(pop_bin,n_parent){
  parents = list()
  p = sample.int(length(pop_bin),16,replace = TRUE)
  for (i in 1:n_parent){
    
  d =abs(max(-f(as.numeric(pop_bin[[p[2*i]]])),-f(as.numeric(pop_bin[[p[2*i-1]]]))) 
    - min(-f(as.numeric(pop_bin[[p[2*i]]])),-f(as.numeric(pop_bin[[p[2*i-1]]]))))/2
  rand = runif(1)*100
  p1=p[c(2*i-1,2*i)]
  if (rand < 50+ d ) { parents[i] = pop_bin[p1[which.max(-f(as.numeric(pop_bin[[p1]])))]] }
  else { parents[[i]] = pop_bin[[p1[which.min(-f(as.numeric(pop_bin[[p1]])))]]] }
  }
  parents
}

###mutation & new generation
mutate = function(parents){ #cross over
  id = sample.int(length(parents),length(parents))
  for ( j in 1:(length(pop_bin)-1)){
    parent1 = parents[id[(2*j)-1]]
    parent2 = parents[id[2*j]]
    c = sample.int(6,1)
    child1=as.numeric(c(parent1[[1]][1:c],parent2[[1]][-(1:c)]))
    child2= as.numeric(c(parent2[[1]][1:c],parent1[[1]][-(1:c)]))
    ch1=child1[1]*(2^5)+child1[2]*(2^4)+child1[3]*(2^3)+child1[4]*(2^2)+child1[5]*2+child1[6] 
    ch2=child2[1]*(2^5)+child2[2]*(2^4)+child2[3]*(2^3)+child2[4]*(2^2)+child2[5]*2+child2[6] 
    
    if (f(ch1) <= f(ch2)) 
      {next.gen[[j+1]] <- child1} 
    else 
      {next.gen[[j+1]] <- child2}
  }
  return(next.gen)
}

################### 
#GENETIC ALGORITHM#
###################
Xs= seq(from = 1, to = 64, length.out = 100)
Ys = f(Xs)
plot(Xs,Ys, type = "l")
points(x=populations[1,],y=f(populations[1,]), pch=20, col="blue")
for(k in 2:M){
  next.gen = list() #retain fittest chromosome 
  next.gen[1]= pop_bin[which.max(fitness[k-1,])]
  parents = selection(pop_bin,2*(length(pop_bin)-1))
  next.gen = mutate(parents)
  populations[k,] = bintonum(next.gen) # save new generation 
  fitness[k,]=-f(populations[k,]) # save fitness vals
  points(x=populations[k,],y=f(populations[k,]), pch=20, col="red")
  pop_bin= next.gen
  readline("next generation")
}
best_answer = populations[M,][which.max(fitness[k-1,])]
cat("Best Answer=", best_answer)
#######################3
####################3
a_baby = (a_d+a_m)/2 + rtrunc(1, spec="cauchy", 
                              a=(30-a_d-a_m)/(2*(a_d-a_m)) ,b=(50-a_d-a_m)/(2*(a_d-a_m)))*(a_d-a_m)/2 
b_baby = (b_d+b_m)/2 + rtrunc(1, spec="cauchy",
                              a=(-b_d-b_m)/(2*(b_d-b_m)) ,b=(10-b_d-b_m)/(2*(b_d-b_m)))*(b_d-b_m)/2 
g_baby = (g_d+g_m)/2 + rtrunc(1, spec="cauchy",
                              a=(-g_d-g_m)/(2*(g_d-g_m)) ,b=(2-g_d-g_m)/(2*(g_d-g_m)))*(g_d-g_m)/2 

a_d = alphas[which.max(alphas)]
a_m = alphas[which.min(alphas)]
b_d = betas[which.max(betas)]
b_m = betas[which.min(betas)]
g_d = gammas[which.max(gammas)]
g_m = gammas[which.min(gammas)]



#################################
C =  rtrunc(1, spec="cauchy", a=-1 ,b=1)
a_baby = (alphas[1]+alphas[2])/2 + C*(alphas[1]-alphas[2])/2 
b_baby = (betas[1]+betas[2])/2 + C*(betas[1]-betas[2])/2 
g_baby = (gammas[1]+gammas[2])/2 + C*(gammas[1]-gammas[2])/2 