#函數
BEMM <- function(rankings, initial.method="mean", it.max=20){
  
  
  #######################################
  #
  #step 1. input the data
  #
  #######################################
  
  
  
  
  #######################################
  #
  #step 2. initialize the parameters
  #
  #######################################
  
  num.ranker = ncol(rankings)
  num.item   = nrow(rankings)
  
  
  #2.1 initialize pi0
  
  score=rep(0,num.item)
  if(initial.method=="mean"){score=apply(rankings, 1, mean)}
  if(initial.method=="median"){
    for(i in 1:num.item){score[i]=quantile(rankings[i,],0.5)}
  }        
  if(initial.method=="geometric"){
    for(i in 1:num.item){score[i]=(prod(rankings[i,]))^(1/num.ranker)}
  }   
  if(initial.method=="random"){
    score=sample(1:num.item,num.item,replace=F)
  } 
  
  
  op.pi0=rank(score,ties.method="random")
  
  #2.2 initialize phi, alpha, omega and phi.h1, tau
  
  op.phi=0.4
  op.alpha=rep(0.4,num.ranker)
  op.omega=rep(0.4,num.ranker)
  
  ########################################################################
  # 
  # Step 3. find the Bayesian estimates of the parameters using Gibbs sampling
  #
  ########################################################################
  
  #to store mcmc samples
  pi0_sample <- matrix(NA, num.item, it.max)
  phi_sample <- c()
  alpha_sample <- matrix(NA, num.ranker, it.max)
  omega_sample <- matrix(NA, num.ranker, it.max)
  
  #3.1 define the likelihood function
  EMM.likelihood <- function(data, tau0, phi, alpha, omega){
    
    #檢查參數是否超出定義域
    domain = (phi>1 |phi<0)
    for(al in alpha){
      domain = (domain |(al>1 | al<0))
    }
    for(om in omega){
      domain = (domain |(om>1 | om<0))
    }
    
    if(domain){
      return(0)
    }else{
      #沒有超出的話才開始計算概似函數
      prob = 1
      
      for(k in 1:ncol(data)){
        tau0_try = tau0
        for(i in 1:nrow(data)){
          index = which(data[,k]==i)
          v = tau0_try[index]
          phi_i = phi*(1-alpha[k]^i)
          part1 = omega[k]*(phi_i^(v-1))/sum(phi_i^(0:(nrow(data)-i)))
          part2 = (1-omega[k])/(nrow(data)-i+1)
          part = part1 +part2
          prob = prob * part
          tau0_try[tau0_try > v] = tau0_try[tau0_try > v]-1
        }
      }
      
      return(prob)
    }
    
  }
  
  #3.2 GIBBS Sampling
  for(it in 1:it.max){
    #3.2.1 update pi0
    
    #proposal為任意交換兩個entities的排名
    op.pi0.pro = op.pi0
    index = sample(1:num.item,2,replace = F)
    a = op.pi0.pro[index[1]]
    b = op.pi0.pro[index[2]]
    op.pi0.pro[index[1]] = b
    op.pi0.pro[index[2]] = a
    
    a = runif(1)
    prob = EMM.likelihood(data=rankings, tau0=op.pi0.pro, phi=op.phi, alpha=op.alpha, omega=op.omega)
    prob = prob/EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha, omega=op.omega)
    prob = min(1,prob)
    if(a<prob){
      op.pi0 = op.pi0.pro
    }else{
      op.pi0 = op.pi0
    }
    
    pi0_sample[,it] <- op.pi0 #紀錄本次的mcmc sample
    
    #3.2.2 update phi
    op.phi.pro = rnorm(1,op.phi,0.05)
    a = runif(1)
    prob = EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi.pro, alpha=op.alpha, omega=op.omega)
    prob = prob/EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha, omega=op.omega)
    prob = min(1,prob)
    if(a<prob){
      op.phi = op.phi.pro
    }else{
      op.phi = op.phi
    }
    
    phi_sample <- c(phi_sample,op.phi)
    
    #3.2.3 update alpha
    op.alpha.pro = op.alpha
    for(i in 1:num.ranker){
      op.alpha.pro[i] = rnorm(1,op.alpha[i],0.05)
      a = runif(1)
      prob = EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha.pro, omega=op.omega)
      prob = prob/EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha, omega=op.omega)
      prob = min(1,prob)
      if(a<prob){
        op.alpha = op.alpha.pro
      }else{
        op.alpha = op.alpha
      }
    }
    
    alpha_sample[,it] <- op.alpha
    
    
    #3.2.4 update omega
    op.omega.pro = op.omega
    for(i in 1:num.ranker){
      op.omega.pro[i] = rnorm(1,op.omega[i],0.05)
      a = runif(1)
      prob = EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha, omega=op.omega.pro)
      prob = prob/EMM.likelihood(data=rankings, tau0=op.pi0, phi=op.phi, alpha=op.alpha, omega=op.omega)
      prob = min(1,prob)
      if(a<prob){
        op.omega = op.omega.pro
      }else{
        op.omega = op.omega
      }
    }
    
    omega_sample[,it] <- op.omega
    
  }
  
  
  return(list(phi_sample=phi_sample,
              omega_sample=omega_sample, 
              alpha_sample=alpha_sample, 
              pi0_sample=pi0_sample))
  
}




#測試資料集
set.seed(2024)
rankings = t(rmm(n=20 , 1:10, theta=1, dist.name = "kendall"))
rankings = rEMM(n=6, tau0=1:30, phi=0.1, alpha=rep(0.5,6), omega=rep(0.8,6))

#測試
set.seed(2024)
out = BEMM(rankings=rankings, initial.method="random", it.max=6000)

#samples of phi
plot(1:6000, out$phi_sample, type = 'l',ylim = c(0,1),xlab = 'iteration',ylab = 'phi')
hist(out$phi_sample[3001:6000],xlim = c(0,1),main = 'phi',xlab = 'MCMC samples of phi',breaks = 20)
(phi.estimate = mean(out$phi_sample[3001:6000]))
abline(v=phi.estimate,col='red')
abline(v=0.1,col='blue')
legend(0.7, 350, legend = c("true phi", "phi estimate"),
       col = c("blue", "red"), lty = 1:1)

#alpha1
plot(1:6000, out$alpha_sample[1,], type = 'l',ylim = c(0,1),xlab = 'alpha1',ylab = 'iteration')
hist(out$alpha_sample[1,3001:6000])
(mean(out$alpha_sample[1,3001:6000]))

#omega
plot(1:6000, out$omega_sample[1,], type = 'l',ylim = c(0,1),xlab = 'omega1',ylab = 'iteration')
hist(out$omega_sample[1,3001:6000])
(mean(out$omega_sample[1,3001:6000]))

#samples of ranker(pi0)
library(tidyverse)
pi0 = out$pi0_sample[,3001:6000]
(agg.ranker = rank(apply(pi0,1,mean)))  #aggregated ranker

pi0.new.rank=pi0.new.name = c() #繪圖
for(i in 1:nrow(pi0)){
  pi0.new.rank = c(pi0.new.rank , pi0[i,])
  pi0.new.name = c(pi0.new.name , rep(paste(as.character(i)), ncol(pi0)))
}
pi0.new = data.frame(
  rank = pi0.new.rank,
  entity = pi0.new.name
)
ggplot(data=pi0.new)+geom_boxplot(aes(x=entity,y=rank))


