# Network Modeling - HS 2021
# C.Stadtfeld, A. Uzaheta, K. Mepham, V.Amati 
# Assignment 3 - Task 1

# Team: Huang Yingshan, Kahlbacher Fabian, Wattenberg Yannick, Aggeler Samuel

#*******************************************************************************
# Task 1.1                                                                  ---- 
# The function "simulation" simulates the network evolution between 
# two time points. 
# Given the network at time t1, denoted by x1, the function simulates the 
# steps of the continuous-time Markov chain defined by a SAOM with outdegree and 
# reciprocity statistics. Unconditional simulation is used.
# The function returns the network at time t2, denoted by x2.
# The structure of the algorithm is describe in an appendix of the slides of
# practical 3.
#*******************************************************************************

#Helper functions to compute P(i->j;x,beta) for a given network and parameters

# Compute outdegree for a given node i in the network x
compute_outdegree <- function(i, x) {
  return(sum(x[i,]))
}

# Compute reciprocity for a given node i in the network x
compute_reciprocity <- function(i, x) {
  sum <- 0
  nvertices <- nrow(x)
  for (j in 1:nvertices) {
    sum <- sum + x[i,j]*x[j,i]   
  }
  return(sum)
}

# Compute the objective function f for a given node i in the network x and 
# weights beta1 and beta2

objective_function <- function(i, x, beta1, beta2) {
  return(beta1 * compute_outdegree(i,x) +  beta2 * compute_reciprocity(i,x))
}

probability_change <- function(i, x, beta1, beta2) {
  n <- nrow(x)
  p <- rep(0., n)
  sum <- 0
  for (j in 1:n) {
    if(i != j) 
        x[i,j] <- !x[i,j] # swap tie i -> j
    p[j] <- exp(objective_function(i, x, beta1, beta2))
    sum <- sum + p[j]
    x[i, j] <- !x[i,j] # restore tie i -> j to original state
  }
  return(p/sum)
}
#' Simulate the network evolution between two time points
#'
#' @param n number of actors in the network
#' @param x1 network at time t1
#' @param lambda rate parameter
#' @param beta1 outdegree parameter
#' @param beta2 reciprocity parameter
#'
#' @return network at time t2
#'
#' @examples
#' netT1 <- matrix(
#'   c(0, 1, 0, 0, 0,
#'     0, 0, 0, 1, 0,
#'     0, 0, 0, 1, 1,
#'     1, 0, 1, 0, 0), 
#'   nrow = 5, ncol =  5, byRow = TRUE)
#' netT2 <- simulation(5, netT1, 4, -2, 0.5)
simulation <- function(n, x1, lambda, beta1, beta2) {
  t <- 0 # time
  x <- x1
  while (t < 1) {
    dt <- rexp(1, n * lambda)
    i <- sample(x=1:n, size=1)
    j <- sample(x=1:n, size=1, prob=probability_change(i, x, beta1, beta2))
    
    if(i != j){
      x[i, j] <- !x[i,j] #flip chosen tie
    }
    # if j == i do nothing
    t = t + dt
  }
  return(x) # return x2
}

#*******************************************************************************
# Task 1.2                                                                  ----
# Consider the two adjacency matrices in the files net1.csv and net2.csv.
# Estimate the parameters of the SAOM with outdegree and reciprocity statistics
#  using the function `siena07`
#*******************************************************************************
library(RSiena)

net1 <- as.matrix(read.csv('net1.csv', header=FALSE))
net2 <- as.matrix(read.csv('net2.csv', header=FALSE))

friendship <- sienaDependent(array(c(net1, net2), dim = c(22, 22, 2)))

mydata <- sienaDataCreate(friendship)
myeff <- getEffects(mydata)

# Specifying the parameter of the algorithm
myAlgorithm <- sienaAlgorithmCreate(
  projname = "net1_net2",
  nsub = 4, n3 = 3000, seed = 1908
)

model0 <- siena07(myAlgorithm,
                  data = mydata, effects = myeff,
                  returnDeps = TRUE, useCluster = TRUE, nbrNodes = 3, batch = FALSE
)

model0

#                             Estimate   Standard   Convergence 
#                                         Error      t-ratio   
#Rate parameters: 
#  0       Rate parameter       4.4326  ( 0.7784   )             
#
#Other parameters: 
# 1. eval outdegree (density)  -1.3178  ( 0.1675   )   0.0436    
# 2. eval reciprocity           1.7292  ( 0.2840   )   0.0156    
#
#Overall maximum convergence ratio:    0.0518 





#*******************************************************************************
# Task 1.3                                                                  ----
# Conditioning on the first observation, generate 1,000 simulations of the 
# network evolution
# Compute the cumulative indegree and outdegree distribution for each 
# simulated network.
# Save the results in a matrix, named indegDist/outDist, in which rows are
# the simulations and columns are the type of triads.
#*******************************************************************************

#' Helper function: Computes the cumulative degree distribution 
#'
#' @param network input adjacency matrix 
#' @param type string indicating the type of degree distribution being computed.
#'   Available types: "indegree" or "outdegree" 
#' @param maxDeg integer with the highest degree value consider
#'
#' @return a integer vector with the cumulative degree distribution
#'
#' @examples
#' net <- matrix(
#'   c(0, 1, 0, 0, 0,
#'     0, 0, 0, 1, 0,
#'     0, 0, 0, 1, 1,
#'     1, 0, 1, 0, 0), 
#'   nrow = 5, ncol =  5, byRow = TRUE)
#' degreeDistribution(net, type = "indegree")
#' degreeDistribution(net, type = "outdegree")
degreeDistribution <- function(network, type = c("indegree", "outdegree"),
                               maxDeg = 8) {
  if (type == "indegree") {
    nodeDegree <- colSums(network)
  } else {
    nodeDegree <- rowSums(network)
  }
  
  possibleValues <- seq(0, min(ncol(network) - 1, maxDeg))
  degreeDist <- sapply(possibleValues, function(x) sum(nodeDegree == x))
  cumulativeDegree <- cumsum(degreeDist)
  names(cumulativeDegree) <- sprintf('Deg%02d', possibleValues)
  return(cumulativeDegree)
}


eval_outdegree <- -1.3178
eval_reciprocity <- 1.7292

N <- 1000

simulate_get_dist <- function(run){
  if(run %% 50 == 0)
    print(run)
  x2 <- simulation(
    n=nrow(net1),
    x1=net1,
    lambda=model0$rate,
    beta1=eval_outdegree,
    beta2=eval_reciprocity
  )
  
  indegDist <- degreeDistribution(
    network=x2,
    type = "indegree",
  )
  
  outdegDist <- degreeDistribution(
    network=x2,
    type = "outdegree",
  )
  
  return(matrix(c(indegDist, outdegDist), nrow = 2, ncol =  9))
}
res <- lapply(seq(1,N), simulate_get_dist)

indegDist <- matrix(NA, 1000, 9)
outdegDist <- matrix(NA, 1000, 9)
for(i in 1:N){
  indegDist[i,] <- res[[i]][1,]
  outdegDist[i,] <-res[[i]][2,]
}

just_for_names <- degreeDistribution(
  network=x2,
  type = "outdegree",
  maxDeg = 8
)

colnames(indegDist) <- names(just_for_names)
colnames(outdegDist) <- names(just_for_names)


#*******************************************************************************
# Task 1.4                                                                  ----
# Fill out the missing part and run the code to obtain the violin plots
#*******************************************************************************
# install.packages(c("tidyverse", "ggplot2"))  # # run this line to install 
library(tidyverse)
library(ggplot2)

# Given the array indegDist, create a data frame from it in a long format, 
# do the same for the observed network statistics at time t2.
# Named the data frame "indegDistDf" and "obsIndegData"

indegDistDf <- data.frame(indegDist) %>% 
  select(where(~ var(.) > 0)) %>%  # Drop statistics with zero variance
  pivot_longer(
    starts_with("Deg"),
    names_to = "degree", names_pattern = "Deg(.+)",
    values_to = "nnodes"
  ) 

# Compute the statistics of the observed network at time t2
obsIndeg <- degreeDistribution(
  network=net2,
  type = "indegree",
  maxDeg = 8
)
  
obsIndegData <- data.frame(
  degree = str_extract(names(obsIndeg), "\\d+"),
  nnodes = obsIndeg) %>% 
  filter(degree %in% unique(indegDistDf$degree)) 

# The following code computes the 5% and the 95% quantiles
# of the distribution of the number of nodes by degree
percIndeg <- indegDistDf %>% 
  group_by(degree) %>% 
  summarise(
    quant05 = quantile(nnodes, prob = 0.05),
    quant95 = quantile(nnodes, prob = 0.95)) %>% 
  pivot_longer(
    starts_with("quant"),
    names_to = "quant", names_pattern = "quant(.+)",
    values_to = "nnodes")


# The following code produces the violin plots
ggplot(indegDistDf, aes(degree, nnodes)) +
  geom_violin(trim = FALSE, scale = "width") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_boxplot(width = 0.2, fill = "gray", outlier.shape = 4) +
  geom_point(data = obsIndegData, 
             col = "red", size = 2) +
  geom_line(data = obsIndegData, aes(group = 1),
            col = "red", size = 0.5) +
  geom_line(data = percIndeg, mapping = aes(group = quant),
            col = "gray", linetype = "dashed") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_x_discrete(labels = as.numeric)


# Mahalanobis and p-value computation                                       ----
# Compute the variance-covariance matrix from the statistics of the 
# simulated networks.
# Compute the inverse of this variance-covariance matrix.
# Mean-center the statistics of the simulated networks.
# Center the statistics of the observed network at time t2 using the means 
# of the statistics of the simulated networks.
# Compute the Mahalanobis distance using the mhd function for 
# the statistics of the simulated networks and the observed network.
# Compute the proportion of simulated networks where the distance is 
# equal or greater than the distance in the observed network.

#' Compute the Mahalanobis distance
#'
#' @param centeredValue numerical vector with the mean centered values
#' @param invCov numerical matrix with the inverse of the variance-covariance
#'   matrix
#'
#' @return numeric value with the Mahalanobis distance of centered value
#'
#' @examples
#' mhd(c(2, 4) - c(1.5, 2), solve(matrix(c(1, 0.8, 0.8, 1), ncol = 2)))
mhd <- function(centeredValue, invCov) {
  t(centeredValue) %*% invCov %*% centeredValue
}

# ---MISSING---

inv_cov <- solve(a=cov(indegDist))

degree_means <- apply(X=indegDist, MARGIN=2, mean)
mhd_indegree_sim <- rep(NA, N)
for (i in 1:N) { #i=1
  mhd_indegree_sim[i] <- mhd(
    centeredValue=(indegDist[i,] - degree_means),
    invCov=inv_cov
  )
}

centeredValues_obs <- 
mhd_indegree_obs <- mhd(
  centeredValue=(obsIndeg - degree_means),
  invCov=inv_cov
)

indeg_precentage <- length(mhd_indegree_sim[mhd_indegree_sim  >= mhd_indegree_obs[1,1]])/N

# So far, we have created the violinplot and the test on the 
# Mahalanobis distance for the auxiliary statistic indegree.
# Now, create the code to generate the violin plot and the test on the 
# Mahalanobis distance for the auxiliary statistic outdegree

# ---MISSING---

outdegDistDf <- data.frame(outdegDist) %>% 
  select(where(~ var(.) > 0)) %>%  # Drop statistics with zero variance
  pivot_longer(
    starts_with("Deg"),
    names_to = "degree", names_pattern = "Deg(.+)",
    values_to = "nnodes"
  ) 

# Compute the statistics of the observed network at time t2
obsOutdeg <- degreeDistribution(
  network=net2,
  type = "outdegree",
  maxDeg = 8
)

obsOutdegData <- data.frame(
  degree = str_extract(names(obsOutdeg), "\\d+"),
  nnodes = obsOutdeg) %>% 
  filter(degree %in% unique(OutdegDistDf$degree)) 

# The following code computes the 5% and the 95% quantiles
# of the distribution of the number of nodes by degree
percOutdeg <- outdegDistDf %>% 
  group_by(degree) %>% 
  summarise(
    quant05 = quantile(nnodes, prob = 0.05),
    quant95 = quantile(nnodes, prob = 0.95)) %>% 
  pivot_longer(
    starts_with("quant"),
    names_to = "quant", names_pattern = "quant(.+)",
    values_to = "nnodes")


# The following code produces the violin plots
ggplot(outdegDistDf, aes(degree, nnodes)) +
  geom_violin(trim = FALSE, scale = "width") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_boxplot(width = 0.2, fill = "gray", outlier.shape = 4) +
  geom_point(data = obsOutdegData, 
             col = "red", size = 2) +
  geom_line(data = obsOutdegData, aes(group = 1),
            col = "red", size = 0.5) +
  geom_line(data = percOutdeg, mapping = aes(group = quant),
            col = "gray", linetype = "dashed") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_x_discrete(labels = as.numeric)


inv_cov <- solve(a=cov(outdegDist))
out_degree_means <- apply(X=outdegDist, MARGIN=2, mean)
mhd_outdegree_sim <- rep(NA, N)
for (i in 1:N) { #i=1
  mhd_outdegree_sim[i] <- mhd(
    centeredValue=(outdegDist[i,] - out_degree_means),
    invCov=inv_cov
  )
}

mhd_outdegree_obs <- mhd(
  centeredValue=(obsOutdeg - degree_means),
  invCov=inv_cov
)

outdeg_precentage <- length(mhd_outdegree_sim[mhd_indegree_sim  >= mhd_outdegree_obs[1,1]])/N
indeg_precentage
outdeg_precentage












































