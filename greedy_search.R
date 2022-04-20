#### 0. Libraries, functions ####
rm(list=ls())
source("netowork-functions.R")

library(igraph)
library(dplyr)
library(zoo) # locf function
library(tidyverse)
library(pbapply)

#### 1. Algorithm function ####

accumulate_seeds = function(x, best.seeds) {
  return(c(x, best.seeds))
}

greedy_search <- function(n_best, model="threshold", network.neighbors, 
                          network.size, ic_iter=10, print_res=F) {
  # Store the best seeds in a vector
  best.seeds <- c()
  
  # Run each iteration on a set of potential seeds. We exclude the best seeds
  # from potential seeds list after each iteration.
  potential_seeds <- c(1:net.size)
  
  test_seeds <- potential_seeds
  
  # Define the number of best seeds you would like to search for.
  if (model == "threshold") {
    
    for (n in 1:n_best) {
      
      # Store the best seed node for each iteration for the given model
      best.seeds[n] <- potential_seeds[which.max(sapply(
        test_seeds, threshold_model, network.neighbors, network.size
      ))]
      
      if (print_res){print(best.seeds[[n]])}
      
      # Subtract best seeds from potential seeds for the next iteration.
      potential_seeds <- potential_seeds[! potential_seeds %in% best.seeds[[n]]]
      
      # To each element of 
      test_seeds = lapply(potential_seeds, accumulate_seeds, best.seeds)
      
    }
  } else if (model == "ic") {
    
    for (n in 1:n_best) {
      
      # Store the best seed node for each iteration for the given model.
      best.seeds[n] <- potential_seeds[which.min(pbapply::pbsapply(
        test_seeds,
        function(node_seed, ic_iter, net.neighbors, net.size) {
          mean(replicate(ic_iter, IC_2(node_seed, net.neighbors, net.size)))
        },
        ic_iter, network.neighbors, network.size
      ))]
      
      if (print_res){print(best.seeds[[n]])}
      
      # Subtract best seeds from potential seeds for the next iteration.
      potential_seeds <- potential_seeds[! potential_seeds %in% best.seeds[[n]]]
      test_seeds = lapply(potential_seeds, accumulate_seeds, best.seeds)
    }
  }
  
  return(best.seeds)
}

#### 2. Run the algorithm on Threshold and IC models ####

net <- readRDS("classnetwork_2022.rds")
net.attributes <- readRDS("class2022_attributes.rds")
net.size = igraph::gorder(net)
thres_attr <- "tRecipe"
net.neighbors <- get_neighbors_list(net)

greedy_search(5, "threshold", net.neighbors, net.size, 
              ic_iter = 10, print_res=T)

greedy_search(5, "ic", net.neighbors, net.size, 
              ic_iter = 100, print_res=T)


