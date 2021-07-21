##################################
## Tools for chapter 4 ##
##################################

# Simulation study to compare UPCS, ACS, SCPS, PoSA, AWS, CPoSA #

#setwd('C:/Users/Emanuela/Dropbox/Furfaro_PhDThesis/Code/Thesis/Chapter4')

## Needed libraries ##
if (!require("spatstat")) install.packages("spatstat")
#library(spatstat) # needed for runifdisc/ppp objects/quadratcount
#library(BalancedSampling) # needed for scps
if (!require("BalancedSampling")) install.packages("BalancedSampling")
library(plyr)
library(dplyr) 
#library(data.table)
if (!require("data.table")) install.packages("data.table")
#library(intergraph)
if (!require("intergraph")) install.packages("intergraph")
#library(magrittr)
if (!require("magrittr")) install.packages("magrittr")
#library(igraph)
if (!require("igraph")) install.packages("igraph")
#library(network)
if (!require("network")) install.packages("network")
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
# developed with Federico Andreis

# UPCS

## pareto / CP sampling functions
par.s <- function(x) {
  
  nn <- sum(x)
  uu <- runif(length(x))
  qq <- uu/(1-uu)/(x/(1-x))
  s <- order(qq)[1:nn]
  
  tmp.s <- numeric(length(x))
  tmp.s[s] <- 1
  
  d <- sum(x*(1-x))
  gma <- (d+x*(1-x))^(-.5)
  g <- (1-x)*sqrt(2*pi)*gma*exp(gma^2*x^2/2)
  B <- sum(sort(g)[1:nn])  
  
  list(cond=(B*sum(tmp.s*g))^-1,s=tmp.s)
  
}

cps.par <- function(x,pareto.s=FALSE) {
  
  tmp <- par.s(x)
  
  s <- tmp$s
  cond <- tmp$cond
  
  if (pareto.s==FALSE) {
    while(runif(1)>cond) {
      
      tmp <- par.s(x)
      s <- tmp$s
      cond <- tmp$cond
      
    }
  }
  
  s
  
}

sample_foo <- function(a) cps.par(a,pareto.s=TRUE) # Pareto Sampling

# ACS
## assign network membership
# From Kisten Sauby code

assignNetworkMembership2 <- function(dataframe, plot.size=1) {
  NetworkID <- x <- NULL
  D <- as.matrix(dist(cbind(dataframe$x, dataframe$y), method="euclidian"))
  D = ifelse(D > plot.size, 0, D)
  D %<>% as.data.frame
  G <- asIgraph(network(D, directed=FALSE))
  dataframe$NetworkID <- clusters(G)$membership
  #-- not Sauby's
  dataframe <- arrange(dataframe, NetworkID)
  network_size <- table(dataframe$NetworkID)
  m <- rep(network_size, network_size)
  dataframe <- as.data.frame(cbind(dataframe,m))
  # dataframe %<>%
  #  group_by(NetworkID) %>%
  #  mutate(m = length(NetworkID)) %>%
  #as.data.frame
  return(dataframe)
}

# partial inclusion probabilities
#pi_i <- function(N, n1, m) {
#  sapply(m, function(m) 
#    1 - exp(
#      sum(log({N - m - n1 + 1} : {N - m})) - sum(log({N - n1 + 1} : N))
#    )
#  )
#}

## Sauby's for creating initial SRS in ACS

createSRSWOR <- function(population, seed, n1) {
  set.seed(seed)
  sample = population[sample(x=1:dim(population)[1], size=n1, replace=F), ]
  # S = merge(population, sample, all.y=TRUE) 	
  sample$Sampling <- "SRSWOR"
  return(sample)		
}

### Sauby's for creating actual ACS

createACS <- function(population, seed=1, n1, y_variable, condition=0, initial_sample=NA) {
  . <- Sampling <- y_val <- NULL
  if (is.data.frame(initial_sample)) {
    S = merge(population, initial_sample, all.y=TRUE) 	
    S$Sampling <- "SRSWOR"
  } else {
    S <- createSRSWOR(population, seed, n1)
  }
  # add the rest of the units for each network in the initial sample
  Z = population %>%
    filter(population$NetworkID %in% S$NetworkID) %>%
    merge(S, all.x=T)
  Networks = filter(Z, eval(parse(text=paste("Z$", y_variable, sep=""))) > condition)
  # if there are units that satisfy the condition, fill in edge units
  if (dim(Networks)[1] > 0) {
    names(Z)[names(Z) == y_variable] <- 'y_val'
    #Z %<>%
    #	as.data.table %>%
    #	setnames(y_variable, "y_val")
    if (dim(Z[which(is.na(Z$Sampling)), ])[1] > 0) {
      Z[which(is.na(Z$Sampling)), ]$Sampling <- "Cluster"
    }
    # fill in edge units
    E = data.table(
      x = as.numeric(rowSums(expand.grid(Networks$x, c(1,-1,0,0)))),
      y = rowSums(expand.grid(Networks$y, c(0,0,1,-1))),
      Sampling = "Edge",
      key = c("x", "y")
    )
    Z %<>% 
      rbind.fill(E) %>% 
      as.data.table %>% 
      setkey("x", "y") %>% 
      unique
    # remove plots outside of population extent
    Z %<>% .[which(Z$x %in% population$x & Z$y %in% population$y)]
    # fill in values for Edge units
    if (dim(Z[ is.na(Z$y_val) ])[1] > 0) {
      Z[ Sampling=="Edge" ]$y_val <- 0 # fill in m
      Z[ y_val==0 & Sampling!="SRSWOR" ]$m <- 0 # fill in m
    }	
    setnames(Z, "y_val", y_variable)
    return(Z)
  } else {
    # if there are NO units that satisfy the condition, stop here and return the SRSWOR sample
    return(Z)
  }
}

# SCPS: no functions needed, in Balanced sampling by A. Grafstrom

# PoSA
# developed with Federico Andreis

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

# AWS

#source('aws_tools.R')

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
### estimators ##
#################
# developed with Federico Andreis

# UPCS
t_HT_UPCS <- function(s, y.s, pik) {
  sum((s*y.s)/pik)
}

# Sauby's - calculate HT estimator for ACS
t_HT_ACS <- function(y, N, n1, pi_i_values=NULL, m=NULL, sampling=NULL, 
                     criterion=NULL) {
  if (!(is.null(sampling)) & !(is.null(criterion))) {
    J = ifelse(y >= criterion | sampling=="SRSWOR", 1, 0)
  } else {
    J = 1
  }
  if (is.null(pi_i_values)) {
    pi_i_values = pi_i(N, n1, m)
  }
  t_HT = sum(y*J/pi_i_values, na.rm=T)/N
  return(t_HT)	
}

# SPCS
t_HT_SPCS <- function(y.s, s.pik) {
  sum(y.s/s.pik)
}

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

# AWS: in tools_aws

# CPoSA

t_cposa <- function(s.weights,s.y) {
  # input:
  # s.weigts = diagonal of the matrix pp (output of cposa)
  # s.y = vector of y for sampled units
  # output: estimate for the total
  
  sum(s.y/s.weights)
  
}

# May 14th

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

# June 25th
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
    # pi_i[3] = pi_init[3]+pi_init[1]*(1-pi_init[2])*d[1]*(1-pi_init[3])*d[2]+
    #                     +pi_init[2]*(1-pi_init[3])*d[2]

# pi_i[4] = pi_init[4]+pi_init[1]*(1-pi_init[2])*d[1]*(1-pi_init[3])*d[2]*d[1]*(1-pi_init[4])*d[3]+
    #                 +pi_init[2]*(1-pi_init[3])*d[2]*d[1]*(1-pi_init[4])*d[3]
    #                 +pi_init[3]*(1-pi_init[4])*d[3]
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

# example
#N = 10
#pi_init = rep(0.2, N)
#y = c(0,1,0,1,0,0,0,1,1,0)
#pi_i = calculate_pi_i(pi_init = pi_init, y=y, cond=0)
#pi_ik = calculate_pi_ik(pi_init = pi_init, y = y, cond = 0, pi_i = pi_i)
# posa_ht_est = numeric(100000)
# posa_pht_est = numeric(100000)
# var_ht_posa = numeric(100000)
# var_pht_posa = numeric(100000)
# for (i in 1:100000) {
#   s = posa(pik = pi_init, Y=y, cond=0)
#   posa_ht_est[i] = t_posa_ht(s.pi_i = pi_i[s==1], s.y = y[s==1])/length(y)
#   posa_pht_est[i] = t_posa_pht(s = s, y.s = y[s==1], pik = pi_init, cond=0)/length(y)
#   var_pht_posa[i] = calculate_var_pht_posa(pi_init=pi_init, pi_i=pi_i, y=y, cond=0, length(y))[[1]]
#   var_ht_posa[i] = calculate_var_ht_posa(pi_i=pi_i, pi_ik = pi_ik, y=y, length(y))
# }
# (mean(posa_ht_est)-sum(y)/length(y))/(sum(y)/length(y))*100
# (mean(posa_pht_est)-sum(y)/length(y))/(sum(y)/length(y))*100
# (mean(unlist(var_pht_posa))-var(posa_pht_est))/(var(posa_pht_est))*100
# (mean(var_ht_posa)-var(posa_ht_est))/(var(posa_ht_est))*100
# var(posa_ht_est)
# var(posa_pht_est)

