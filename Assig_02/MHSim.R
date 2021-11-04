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
  # calculates the number of edges
  return(sum(net))
}

z2 <- function(net){
  # calculates the number of reciprocal dyads
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
  # calculating the terms needed to calculate the probability of the next state 
  statistic <- exp(theta1 * z1(net) + theta2 * z2(net))
  
  net[tie[1], tie[2]] <- !net[tie[1], tie[2]] # flip a tie
  flipped_tie_statistic <- exp(theta1 * z1(net) + theta2 * z2(net))
  # ***************************************************************
  ## ask whether this formula is the one we should use???
  # ***************************************************************
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  #                --- MISSING---
  p_next_state <- flipped_tie_statistic/(statistic + flipped_tie_statistic)
  
  # Select the next state: 

  #                --- MISSING---
  
  # In this section it is evaluated whether the tie should be flipped or not.
  # Since we already flipped the tie we are doing exactly the opposite and 
  # evaluate whether it should be flipped back to the original state
  # (= not flipping at all) or not (=flipping). 
  
  if (runif(1) > p_next_state) {
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
  
  # Perform the burning steps
  #                --- MISSING---
  # iterate over the given amount of burning steps and execute 
  # Metropolis-Hastings steps
  for (burninStep in 1:burnin) {
    net <- MHstep(net, theta1, theta2)
  }

  # After the burning phase we draw the networks
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
  
  # Return the simulated networks and the statistics
  list(netSim=netSim, statSim=statSim)
}


