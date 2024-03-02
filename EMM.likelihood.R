#likelihood of emm

EMM.likelihood <- function(data, tau0, phi=0.3, alpha, omega){
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

#測1
EMM.likelihood(data=rankings,
               tau0 = 1:100,
               alpha=rep(0.3,10),
               omega = rep(0.3,10))


#測試2

data = list(
  matrix(c(1,2,3),ncol = 1),
  matrix(c(1,3,2),ncol = 1),
  matrix(c(2,1,3),ncol = 1),
  matrix(c(2,3,1),ncol = 1),
  matrix(c(3,1,2),ncol = 1),
  matrix(c(3,2,1),ncol = 1)
)

tau0 = c(1,2,3)


data = list(
  matrix(c(2,1),ncol = 1),
  matrix(c(1,2),ncol = 1)
)

tau0 = c(1,2)

total <- 0
for(i in 1:6){
  total <- total+EMM.likelihood(data=data[[i]],
                                tau0 = 1:6,
                                alpha=rep(0.3,30),
                                omega = rep(0.3,30))
  print(EMM.likelihood(data=data[[i]],
                       tau0 = 1:6,
                       alpha=rep(0.3,30),
                       omega = rep(0.3,30)))
}
total