#These codes are refer to R package "ExtMallows" (Han Li et al.),for more details, please cheak out 'https://cran.r-project.org/web/packages/ExtMallows/index.html'
BEMM <- function(rankings, initial.method="mean", it.max=20){
  
  ###input:
  #ranking: a matrix, element at row i and col j stands for entity i's rank in ranker j
  #initial.method:initial value of true rank(tau) for mcmc iterations
  #it.max:number of mcmc iterations
  ###output: 
  #phi_sample:mcmc samples of phi
  #pi0_sample:mcmc samples of true ranker
  #alpha_sample:mcmc samples of alpha
  #omega_sample:mcmc samples of omega
  
  
  #######################################
  #
  #step 1. Check the data
  #
  #######################################
  
  num.ranker = ncol(rankings)
  num.item   = nrow(rankings)
  
  
  
  
  #######################################
  #
  #step 2. initialize the parameters
  #
  #######################################
  
  
  
  #2.1 initialize pi0
  
  score=rep(0,num.item)
  
  mean.na <- function(x){
    return(mean(x, na.rm = T))
  }
  
  if(initial.method=="mean"){score=apply(rankings, 1, mean.na)}
    
  if(initial.method=="random"){
    score=sample(1:num.item,num.item,replace=F)
  } 
  
  
  op.pi0=rank(score,ties.method="random")
  
  #2.2 initialize phi, alpha, omega and phi.h1, tau
  
  op.phi=0.5
  op.alpha=rep(0.5,num.ranker)
  op.omega=rep(0.5,num.ranker)
  
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
      
      rank.length = c()
      for(k in 1:ncol(data)){
        rank.length <- c(rank.length, length(data[,k])-sum(is.na(data[,k])))
      }
      
      for(k in 1:ncol(data)){
        tau0_try = tau0
        #i in 1:nrow(data)
        
        for(i in 1:rank.length[k]){
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
rankings = rEMM(n=5, tau0=1:5, phi=0.1, alpha=rep(0.5,5), omega=rep(0.8,5))

#測試
set.seed(2024)
library(PAMA)
NBAFL = cbind(NBANFL()$NBAPL,NBANFL()$NBA)
out = BEMM(rankings=NBAFL, initial.method="mean", it.max=6000)

#samples of phi
plot(1:6000, out$phi_sample, type = 'l',ylim = c(0,1),xlab = 'iteration',ylab = '',main = 'phi')

hist(out$phi_sample[3001:6000],xlim = c(0,1),main = '',xlab = 'MCMC samples of phi',breaks = 20)
(phi.estimate = mean(out$phi_sample[3001:6000]))
abline(v=phi.estimate,col='red')
abline(v=0.1,col='blue')
legend(0.7, 350, legend = c("true phi", "phi estimate"),
       col = c("blue", "red"), lty = 1:1)

#alpha1
bayes.alpha= apply(out$alpha_sample[,3001:6000],1,mean)
alpha.name = c()

par(mfrow=c(2,5))
for(k in 1:nrow(out$alpha_sample)){
  #plot(1:6000, out$alpha_sample[k,], type = 'l',ylim = c(0,1),xlab ='iteration',ylab='',main = paste('alpha',k))
  hist(out$alpha_sample[k,3001:6000],main='',xlab = paste('MCMC samples of alpha',k))
  abline(v = bayes.alpha[k],col='red')
  #alpha.name = c(alpha.name ,paste('alpha',k))
}

out.alpha.new=matrix(bayes.alpha,nrow=1)
colnames(out.alpha.new) = alpha.name
rownames(out.alpha.new) = c('Posterior mean')

out.alpha = c()
alpha.name = c()
lev = c()
for(i in 1:nrow(out$alpha_sample)){
  out.alpha = c(out.alpha, out$alpha_sample[i,3001:6000] )
  alpha.name = c(alpha.name, rep(paste('alpha',as.character(i)),3000))
  lev = c(lev ,paste('alpha',as.character(i)))
}
out.alpha.new = data.frame(
  alpha = out.alpha,
  group = alpha.name
)
out.alpha.new$group = factor(out.alpha.new$group , levels=lev )
boxplot(alpha ~ group, data = out.alpha.new,ylab='')

#omega
bayes.omega = apply(out$omega_sample[,3001:6000],1,mean)
omega.name = c()

par(mfrow=c(2,5))
for(k in 1:nrow(out$omega_sample)){
  plot(1:6000, out$omega_sample[k,], type = 'l',ylim = c(0,1),xlab = 'iteration',main=paste('omega',k),ylab = '')
  #hist(out$omega_sample[k,3001:6000],main='',xlab = paste('MCMC samples of omega',k))
  #abline(v = bayes.omega[k],col='red')
  #omega.name = c(omega.name ,paste('omega',k))
}

out.omega.new=matrix(bayes.omega,nrow=1)
colnames(out.omega.new) = omega.name
rownames(out.omega.new) = c('Posterior mean')

plot(1:6000, out$omega_sample[10,], type = 'l',ylim = c(0,1),xlab = 'omega1',ylab = 'iteration')
hist(out$omega_sample[10,3001:6000])
(mean(out$omega_sample[10,3001:6000]))

out.omega = c()
omega.name = c()
lev = c()
for(i in 1:nrow(out$omega_sample)){
  out.omega = c(out.omega, out$omega_sample[i,3001:6000] )
  omega.name = c(omega.name, rep(paste('omega',as.character(i)),3000))
  lev = c(lev ,paste('omega',as.character(i)))
}
out.omega.new = data.frame(
  omega = out.omega,
  group = omega.name
)
out.omega.new$group = factor(out.omega.new$group , levels=lev )
boxplot(omega ~ group, data = out.omega.new)

#samples of ranker(pi0)
pi0 = out$pi0_sample[,3001:6000]
(agg.ranker = rank(apply(pi0,1,mean)))  #aggregated ranker

pi0.new.rank=pi0.new.name =lev= c() #繪圖

Name = c('Heat','Thunder','Spurs','Celties','Clippers',
         'Lakers','Pacers','76ers','Mavericks','Bulls',
         'Knicks','Grizzlies','Nuggets','Magic','Hawks',
         'Jazz','TrailBlazers','Rockets','Bucks','Suns',
         'Nets','Warriors','Timberwolves','Hornets','Pistons',
         'Kings','Wizards','Raptors','Cavaliers','Bobcats')

for(i in 1:nrow(pi0)){
  pi0.new.rank = c(pi0.new.rank , pi0[i,])
  pi0.new.name = c(pi0.new.name , rep(Name[i], ncol(pi0)))
  lev = c(lev,as.character(i))
}
pi0.new = data.frame(
  rank = pi0.new.rank,
  entity = pi0.new.name
)
pi0.new$entity = factor(pi0.new$entity , levels=Name )
ggplot(data = pi0.new)+
  geom_boxplot(aes(x=entity,y=rank))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

boxplot(rank ~ entity, data = pi0.new)

agg.team = rbind(Name,agg.ranker)
rownames(agg.team) = c('NBA Team','aggregated rank')
