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

MHstep <- function(net, theta1, theta2){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace=FALSE
  tie <- sample(1:nvertices,2,replace=FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  # Compute the change statistics
  
  #                --- MISSING---
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  #                --- MISSING---
  
  # Select the next state: 

  #                --- MISSING---
  
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
  
  # After the burnin phase we draw the networks
  # The simulated networks and statistics are stored in the objects
  # netSim and statSim
  netSim <- array(0,dim=c(nvertices,nvertices,nNet))
  statSim <- matrix(0,nNet,2)
  thinningSteps <- 0 # counter for the number of thinning steps
  netCounter <- 1 # counter for the number of simulated network
  
  #                --- MISSING---
  
  # Return the simulated networks and the statistics
  list(netSim=netSim, statSim=statSim)
}