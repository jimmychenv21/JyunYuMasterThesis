
#此函數可以生成EMM ranking lists
#n:一個正整數，代表要生成幾個排序
#tau0:一個vector，為真實排名
#phi:一個數字，phi
#alpha:一個vector，alpha1~alphak
#omega:一個vector，omega1~omegak

rEMM <- function(n, tau0, phi, alpha, omega){
  
  num.item = length(tau0)
  num.ranker = n
  output = matrix(NA, num.item, num.ranker)
  
  for(k in 1:num.ranker){
    
    tau0.try = tau0
    n0 = num.item
    pass = c()
    
    for(i in 1:num.item){
      
      phi.new = phi*(1-alpha[k]^i)
      prob = phi.new^(1:n0)
      prob = prob/sum(prob)
      
      decide = sample(1:2,1,prob=c(omega[k],1-omega[k]))
      if(decide == 1){ #有omega的機率以EMM決定第i名的entity
        v = sample(1:n0,1,prob = prob )
        index = which(tau0.try == v)
        output[index,k] = i
        
      }else{          #有1-omega的機率隨機決定第i名的entity
        v = sample(1:n0,1)
        index = which(tau0.try == v)
        output[index,k] = i
        
      }
      
      tau0.try[tau0.try==v] = -1
      tau0.try[tau0.try>v] = tau0.try[tau0.try>v]-1
      n0 = n0-1
      
    }
    
    
  }
  
  return(output)
}




#測試

rEMM(n=5, tau0=1:5, phi=0.1, alpha=rep(0.5,5), omega=rep(0.9,5))
