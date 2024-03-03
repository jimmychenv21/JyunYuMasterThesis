#likelihood of EMM

EMM.likelihood <- function(data, tau0, phi, alpha, omega){
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

#測試1(看partial list是否可以運作)
library(PAMA)
NBAFL = NBANFL()$NBAPL

NBAFL=NBAFL[,25:34]
EMM.likelihood(data=NBAFL,
               tau0 = 1:30,
               phi = 0.1,
               alpha=rep(0.8,10),
               omega = rep(0.1,10))


#測試2(機率加總起來是1)

data = list(
  matrix(c(1,2,3),ncol = 1),
  matrix(c(1,3,2),ncol = 1),
  matrix(c(2,1,3),ncol = 1),
  matrix(c(2,3,1),ncol = 1),
  matrix(c(3,1,2),ncol = 1),
  matrix(c(3,2,1),ncol = 1)
)

tau0 = c(1,2,3)


total <- 0
for(i in 1:6){
  total <- total+EMM.likelihood(data=data[[i]],
                                tau0 = tau0,
                                phi=0.1,
                                alpha=rep(0.3,6),
                                omega = rep(0.3,6))
  print(EMM.likelihood(data=data[[i]],
                       tau0 = tau0,
                       phi=0.1,
                       alpha=rep(0.3,6),
                       omega = rep(0.3,6)))
}
total