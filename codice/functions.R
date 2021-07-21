##################################
## Functions ##
##################################
# This code was originally developed with Federico Andreis

# Simulation study to compare PoSA estimators, and CPoSA #

## Needed libraries ##
if (!require("plyr")) install.packages("plyr") 
if (!require("dplyr")) install.packages("dplyr")
if (!require("spatstat")) install.packages("spatstat")
## ---------------- ##

######################################
### For simulating populations #######
# developed with Federico Andreis

nclust2 <- function(x0, y0, radius, kk) {
  n <- rpois(1,kk) # modify to account for between-cluster variability
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

nclust2_ellipse <- function(x0, y0, radius, kk) {
  n <- rpois(1,kk) # modify to account for between-cluster variability
  return(affine(runifdisc(n, radius, centre=c(x0, y0)), mat=diag(c(0.5,1))))
}

############################
##### Sampling designs #####

#### produces PoSA sample with free sample size
posa <- function(pik, Y, cond) {
  n <- sum(pik)
  N <- length(Y)
  Y_new <- (Y>cond)*1
  
  s <- numeric(n)
  
  s[1] <- rbinom(1,1,pik[1])
  
  i <- 2
  
  while(i<=N) {
    
    s[i] <- rbinom(1,1,pik[i])*(1-s[i-1]*Y_new[i-1])+s[i-1]*Y_new[i-1]
    
    i <- i+1
    
  }
  s
}

# CPoSA

cposa <- function(pik, Y, cond) {
  # take a sample with conditional posa design
  # input:
  # pik= vector of inclusion probabilities for the whole population
  # Y= vector of Y in the whole population
  # cond= condition to be greater than!
  # output:
  # the sample pp[N,]=S = (0,0,1,...1, ..,0), and the design weights to be
  # used in the estimation
  
  Y <- (Y>cond)*1
  
  n <- sum(pik)
  N <- length(Y)
  
  s <-  numeric(N)
  
  pp <- matrix(ncol=N, nrow=N)
  
  i <- 0
  
  w <- NULL
  
  w <- 1/(N-1:N+1)
  
  while(i<N) {
    
    i <- i+1
    
    if(i==1) {
      
      pp[i:N,i] <- rbinom(1,1,max(min(pik[i],1),0))
      pp[i,i+1] <- pp[i,i]*Y[i]+(1-pp[i,i]*Y[i])*(pik[i+1]-(pp[i,i]-pik[i])*w[i])
      pp[i,(i+2):N] <- pik[i+1]-(pp[i,i]-pik[i])*w[i]
      
    } else {
      
      if (i<N-1) {
        pp[i:N,i] <- rbinom(1,1,max(0,min(pp[i-1,i],1)))
        pp[i,i+1] <- pp[i,i]*Y[i]+(1-pp[i,i]*Y[i])*(pp[i-1,i+1]-(pp[i,i]-pp[i-1,i])*w[i])
        pp[i,(i+2):N] <- pp[i-1,(i+2):N]-(pp[i,i]-pp[i-1,i])*w[i]
      } else if (i<N) {
        
        pp[i:N,i] <- rbinom(1,1,max(0,min(pp[i-1,i],1)))
        pp[i,i+1] <- pp[i,i]*Y[i]+(1-pp[i,i]*Y[i])*(pp[i-1,i+1]-(pp[i,i]-pp[i-1,i])*w[i])
        
      } else if (i==N) {
        
        pp[i:N,i] <- rbinom(1,1,max(0,min(pp[i-1,i],1)))
        
      }
    }
    
  }
  
  return(list(s=pp[N,], weights=c(pik[1],diag(pp[1:(N-1),2:N]))))
  
}

#################
### Estimators ##
#################

## produces PoSA estimates with free sample size
t_posa <- function(s,y.s,pik,cond) {
  
  NN <- length(pik)
  
  yy <- rep(NA,NN)
  yy[s==1] <- y.s
  yy_new <- rep(NA,NN)
  for (i in 1:NN)
  {
    if (yy[i] >= cond & !is.na(yy[i])) yy_new[i] <- 1 
    if (yy[i] < cond & !is.na(yy[i])) yy_new[i] <- 0
  }
  
  pp <- numeric(NN)
  
  pp[1] <- pik[1]
  for (i in 2:NN) pp[i] <- (s[i-1]==0|(s[i-1]==1&yy_new[i-1]==0))*pik[i]+(s[i-1]==1&yy_new[i-1]==1)*1
  
  pp.s <- pp[s==1]
  
  sum(y.s/pp.s)
  
}

s_weights_posa <- function(s,y.s,pik,cond) {
  
  NN <- length(pik)
  
  yy <- rep(NA,NN)
  yy[s==1] <- y.s
  yy_new <- rep(NA,NN)
  for (i in 1:N)
  {
    if (yy[i] >= cond & !is.na(yy[i])) yy_new[i] <- 1 
    if (yy[i] < cond & !is.na(yy[i])) yy_new[i] <- 0
  }
  
  pp <- numeric(NN)
  
  pp[1] <- pik[1]
  for (i in 2:NN) pp[i] <- (s[i-1]==0|(s[i-1]==1&yy_new[i-1]==0))*pik[i]+(s[i-1]==1&yy_new[i-1]==1)*1
  
  pp.s <- pp[s==1]
  
  return(pp.s)
  
}

# CPoSA

t_cposa <- function(s.weights,s.y) {
  # input:
  # s.weigts = diagonal of the matrix pp (output of cposa)
  # s.y = vector of y for sampled units
  # output: estimate for the total
  
  sum(s.y/s.weights)
  
}

var_posa <- function(s.weights, s.y, N) {
  # input:
  # s.weigts = diagonal of the matrix pp (output of cposa)
  # s.y = vector of y for sampled units
  # N = population size
  # output: estimate for the variance based on formula (11) paper Mecatti, Sismanidis, Furfaro 2020
  
  first_term <- sum(((s.y/s.weights)^2)*(1-s.weights))
  
  # find subsequent values
  ones <- which(s.weights %in% c(1))
  correlated_pairs <- cbind(ones-1, ones)
  
  y_i <- s.y[correlated_pairs[,1]]
  y_i1 <- s.y[correlated_pairs[,2]]
  s.weights_i <- s.weights[correlated_pairs[,1]]
  s.weights_i1 <- s.weights[correlated_pairs[,2]]

  num <- y_i + s.weights_i*(1-y_i)-s.weights_i1
  den <- s.weights_i*s.weights_i1*(y_i+s.weights_i1*(1-y_i))
  
  second_term <- 2*sum(y_i*y_i1*num/den)
  
  var_posa_value <- (1/N^2)*(first_term-second_term)
  
  return(var_posa_value)
}

var_cposa <- function(s.weights, s.y, N) {
  # input:
  # s.weigts = diagonal of the matrix pp (output of cposa)
  # s.y = vector of y for sampled units
  # N = population size
  # output: estimate for the variance based on formula (16) paper Mecatti, Sismanidis, Furfaro 2020
  
  first_term <- sum(((s.y/s.weights)^2)*(1-s.weights))

  # find subsequent values
  ones <- which(s.weights %in% c(1))
  correlated_pairs <- cbind(ones-1, ones)
  
  y_i <- s.y[correlated_pairs[,1]]
  y_i1 <- s.y[correlated_pairs[,2]]
  s.weights_i <- s.weights[correlated_pairs[,1]]
  s.weights_i1 <- s.weights[correlated_pairs[,2]]
  
  num <- 1-s.weights_i1
  den <- s.weights_i*s.weights_i1*(y_i+(1-y_i)*s.weights_i1)
  
  second_term <- 2*sum(y_i^2*y_i1*(num/den))

  var_cposa_value <- (1/N^2)*(first_term - second_term)
  
  return(var_cposa_value)
}

# equation 6, Mecatti-Sismanidis-Furfaro-Conti draft June 2021
calculate_pi_i <- function(pi_init, y, cond) {
  d = ifelse(y>cond, 1, 0)
  M = length(pi_init)
  pi_i = numeric(M)
  pi_i[1] = pi_init[1]
  pi_i[2] = pi_init[2] + pi_init[1]*(1-pi_init[2])*d[1]
  for (i in 3:M) {
    second_term_tmp = 0
    second_term = 0
    for (j in 1:(i-1)) {
      second_term_tmp = pi_init[j]*prod((1-pi_init[(j+1):i])*d[j:(i-1)])
      second_term = second_term + second_term_tmp
    }
    pi_i[i] = pi_init[i]+second_term
  }
  return(pi_i)
}

# Formula 19 with pi_i instead of posa weights
t_posa_ht <- function(s.pi_i,s.y) {
  # input:
  # s.pi_i = vector of pi_i of the sampled units as in equation 6
  # s.y = vector of y for sampled units
  # output: estimate for the total
  sum(s.y/s.pi_i)
}

# at the moment this is the same as t_posa() function above
t_posa_pht <- function(s,y.s,pik,cond) {
  
  M <- length(pik)
  
  yy <- rep(NA,M)
  yy[s==1] <- y.s
  yy_new <- rep(NA,M)
  for (i in 1:M)
  {
    if (yy[i] > cond & !is.na(yy[i])) yy_new[i] <- 1 
    if (yy[i] <= cond & !is.na(yy[i])) yy_new[i] <- 0
  }
  
  pp <- numeric(M)
  
  pp[1] <- pik[1]
  for (i in 2:M) pp[i] <- (s[i-1]==0|(s[i-1]==1&yy_new[i-1]==0))*pik[i]+(s[i-1]==1&yy_new[i-1]==1)*1
  
  pp.s <- pp[s==1]
  
  sum(y.s/pp.s)
  
}

# pi_ik, Eq 8
calculate_pi_ik <- function(pi_init, y, cond, pi_i) {
# returns a matrix
# pi_i: marginal inclusion probabilties
# pi_init: initial inclusion probabilities
# y: pop
# cond: adaptive condition
# M: number of areas, computed within the function
  d = ifelse(y>cond, 1, 0)
  M = length(pi_init)
  pi_ik = matrix(ncol=M,nrow=M)
  for (i in 1:(M-1)) {
    for (k in 1:(M-i)) {
      pi_ik[i,(i+k)] = pi_i[i]*(1-pi_i[i])*prod(1-pi_init[i:(i+k)])*prod(d[i:(i+k-1)]) + pi_i[i]*pi_i[i+k]
    }
  }
  return(pi_ik)
}

# Unnumbered equation after eq. 24
calculate_var_pht_posa <- function(pi_init, pi_i, y, cond, N) {
  # pi_init: initial inclusion probability 
  # pi_i: marginal inclusion probability 
  # y: number of cases in each area
  # cond: adaptive condition 
  # N: number ofunits in the population (!= number of areas)
  d = ifelse(y>cond, 1, 0)
  M = length(pi_init)
  term_1 = sum((1/pi_init-1)*y^2)
  term_2 = sum((1/pi_init[2:M]-1)*d[1:(M-1)]*pi_i[1:(M-1)]*y[2:M]^2)
  var_pht = 1/N^2*(term_1-term_2)
  var_ht_p_init = 1/N^2*term_1
  eff_gain = 1/N^2*term_2
  return(list(var_pht = var_pht, var_ht_p_init = var_ht_p_init, 
              eff_gain = eff_gain))
}


# Equation in AggiuntaSim-toEF-1.pdf
calculate_var_ht_posa <- function(pi_i, pi_init, y, N, cond) {
  # pi_i: marginal inclusion probabilities of all units
  # pi_ik: second order inclusion probabilities of all units (should be a matrix)
  # y: y over all population
  # N: population size
  M = length(pi_i)
  d = ifelse(y>cond, 1, 0)
  term_1 = sum((y^2/pi_i)*(1-pi_i))
  term_2 = 0
  term_2_tmp = 0
  for (i in 1:(M-1)) {
    pi_ik_pi_i_ik = 0
    term_2_2_tmp = 0
    term_2_2 = 0
    for (k in 1:(M-i)) {
    pi_ik_pi_i_ik = pi_i[i]*(1-pi_i[i])*prod(1-pi_init[(i+1):(i+k)])*prod(d[i:(i+k-1)])
    term_2_2_tmp = (y[(i+k)]/pi_i[(i+k)])*(pi_ik_pi_i_ik)
    term_2_2 = term_2_2_tmp + term_2_2
    }
    term_2_tmp = y[i]/pi_i[i]*term_2_2
    term_2 = term_2 + term_2_tmp
  }
  value_var_ht = 1/N^2*(term_1+2*term_2)
  return(value_var_ht)
}

