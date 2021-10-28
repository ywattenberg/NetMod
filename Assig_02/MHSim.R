# Network Modeling - HS 2021
# C. Stadtfeld, V. Amati, A. Uzaheta, K. Mepham  
# Assignment 2 - Task 1

#************************************************************************************
# The function MHstep simulates the the Metropolis-Hastings step that defines
# the Markov chain whose stationary distribution is the ERGM with 
# edge and mutual statistics
#************************************************************************************
# Input:
# net: adjacency matrix of the network
# theta1, theta2: parameters of the ERGM

# Output: next state of the chain

z1 <- function(net){
  return(sum(net))
}

z2 <- function(net){
  sum <- 0
  nvertices <- nrow(net)
  for (i in 1:nvertices) {
    for (j in 1:i) {
      sum <- sum + net[i,j] * net[j,i] 
    }
  }
  return(sum)
}

MHstep <- function(net, theta1, theta2){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace=FALSE
  tie <- sample(1:nvertices,2,replace=FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  # Compute the change statistics
  
  #                --- MISSING---
  stat <- exp(theta1 * z1(net) + theta2 * z2(net))
  
  net[tie[1],tie[2]] <- !net[tie[1],tie[2]] # flip a tie
  change_stat <- exp(theta1 * z1(net) + theta2 * z2(net))
  # ***************************************************************
  ## ask what they exactly want for the change statistics
  # ***************************************************************
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  #                --- MISSING---
  acceptance_ratio <- change_stat/(stat + change_stat)
  
  # Select the next state: 

  #                --- MISSING---
  
  if (runif(1) > acceptance_ratio) {
    net[tie[1], tie[2]] <- !net[tie[1], tie[2]] # flip it back
  }
  # Return the next state of the chain
  return (net)
}



#************************************************************************************
# The function MarkovChain simulates the networks from the ERGM with 
# edge and mutual statistics
#************************************************************************************
# Input:
# net: adjacency matrix of the network
# theta1, theta2: parameters of the ERGM
# burnin: number of steps to reach the stationary distribution
# thinning: number of steps between simulated network 
# nNet: number of simulated networks

# Output:
# 1. netSim: adjacency matrix of the simulated networks
# 2. statSim; value of the edge and mutual statistic

MarkovChain <- function(net, theta1, theta2, burnin=10000, 
                        thinning = 1000, nNet = 1000){
  
  # Burnin phase: repeating the steps of the chain "burnin" times  
  nvertices <- nrow(net)
  burninStep <- 1 # counter for the number of burnin steps
  
  # Perform the burnin steps
  #                --- MISSING---
  for (burninStep in 1:burnin) {
    net <- MHstep(net, theta1, theta2)
  }

  # After the burnin phase we draw the networks
  # The simulated networks and statistics are stored in the objects
  # netSim and statSim
  netSim <- array(0,dim=c(nvertices,nvertices,nNet))
  statSim <- matrix(0,nNet,2)
  thinningSteps <- 0 # counter for the number of thinning steps
  netCounter <- 1 # counter for the number of simulated network
  
  #                --- MISSING---
  for (netCounter in 1:nNet) {
    for(thinningSteps in 1:thinning) {
      net <- MHstep(net, theta1, theta2)
    }
    netSim[,,netCounter] <- net
    statSim[netCounter,] <- c(z1(net), z2(net))
  }
  # ***************************************************************
  ## ask: probabilities or counts??
  # ***************************************************************
  
  # Return the simulated networks and the statistics
  list(netSim=netSim, statSim=statSim)
}

#setwd("C:/Projects/NetMod/Assig_02/")
communication <- as.matrix(read.csv("data/communication.csv", header=F))


net_0 <- matrix(0, 38,38)

sim_net <- MarkovChain(
  net=net_0,
  theta1=-0.6,
  theta2=1.5,
)

sim_net_2 <- MarkovChain(
  net=net_0,
  theta1=-1.737,
  theta2=1.632,
)

mean(sim_net_2$statSim[,1])
mean(sim_net_2$statSim[,2])

mean(sim_net$statSim[,1])
mean(sim_net$statSim[,2])

z1(communication)
z2(communication)



dim(sim_net$netSim)

sim_net$netSim[,,1]

plot(sim_net$statSim[,1], sim_net$statSim[,2], type='l')


# compare with ergm
library(network) 
library(ergm)
netw <- network(communication, directed=TRUE)

ergm(netw ~ edges + mutual) #gwesp(decay=0.3,fixed=TRUE)
#edges    mutual  
#-1.737   1.632