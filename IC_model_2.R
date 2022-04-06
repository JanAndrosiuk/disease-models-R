#### Independent Cascade Algorithm with various probabilities ####

# ASSUMPTION 1:
  # # Case 2 follows differs a bit from the case 1 example.
  # In the previous case, chances of the unactive to become active were 
  # determined by the contagious period of its neighbors and their probability.
  # Here we assume that nodes that became active will stay active till
  # the end of the simulation.
  # Therefore, here unactive node's chances are determined 
  # by the number of node's active neighbors and their probability.

# ASSUMPTION 2:
  # Based on A1: At some time step the network will be 100% infected.

# ASSUMPTION 3:
  # Neighbors' probabilities are not constant.


#### 0. Libraries, functions ####
rm(list=ls())

library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path)) 
# getwd()

source("netowork_functions.R")

# library(visNetwork)
# library(htmlwidgets)
library(igraph)
library(dplyr)
library(zoo) # locf function
library(tidyverse)
library(pbapply)


#### 1. Simulation  function ####

IC_2 <- function(spread.seeds, net.neighbors, net.size) {
  
  # Save number of newly activated nodes for each day in the log.
  spread.log <- c()
  
  # Column of the network attributes containing probabilities of contagion
  # for each node.
  proba_attr <- "pYoutube"
  
  # Day 0
  
  # spread.status informs about nodes that were activated during the simulation.
  # Initially every nodes is inactive (State=0).
  spread.status = rep(0, net.size)
  
  # Day 1
  
  spread.status[spread.seeds] <- 1
  spread.log[[1]] <- length(spread.seeds)
  
  # Day 2 till the end
  
  while (sum(spread.status) != net.size) {
    
    # Get list of active nodes which will have a chance to spread phenomena
    # across their neighbors.
    spread.active = which(spread.status == 1)
    
    # For nodes that are unactive, but have at least one active neighbor, find
    # other active neighbors. Firstly find unactive neighbors of active nodes.
    # Then find active neighbors of those unactive nodes.
    spread.neighbors <- setdiff_nongreedy(
      unlist(net.neighbors[spread.active], use.names=F),
      which(spread.status == 1)
    ) %>%
      net.neighbors[.] |>
      lapply(
        setdiff, which(spread.status == 0)
      )
    
    # Get probabilities of spread "from" each of analyzed neighbor,
    # as these are neighbors that spread contagion to unactive nodes.
    spread.neighbor_probas <- net.attributes[,proba_attr][
      unlist(spread.neighbors, use.names = F)
    ]
    
    # Once we have probabilities for neighbors of unactive nodes,
    # we need to assign those probabilities to each unactive node.
    # That is why I create a vector of repeated unactive nodes, where number of
    # repetitions refers to the number of each unactive node's neighbors. 
    spread.unactive_rep <- rep(
      as.numeric(names(spread.neighbors)),
      sapply(spread.neighbors, length)
    )
    
    # Choose neighbors that will be activated in this iteration.
    spread.new <- unique(spread.unactive_rep[
      runif(length(spread.neighbor_probas), 0, 1) < spread.neighbor_probas
    ])
    
    # If there are any new infected cases, change their status, and save their
    # occurrence to the log.
    if (length(spread.new) > 0) {
      spread.status[spread.new] <- 1
      spread.log[[length(spread.log)+1]] <- length(spread.new)  
    } else {
      spread.log[[length(spread.log)+1]] <- 0
    }
  }
  return(cumsum(spread.log) / net.size)
}


IC_2_average <- function(net.size, no_reps, network_neighbors) {
  return(
    pbapply::pbreplicate(
      no_reps,
      lapply(
        1:net.size, IC_2, network_neighbors, net.size
      )
    )
  )
}


#### 2. Simulation parameters ####

# Load network object
net <- readRDS("data/net1.rds")
net.attributes <- readRDS("data/net1_attr.rds")
net.size = igraph::gorder(net)

# Define start seeds
spread.seeds = c(4)

# nested list of neighborhood for each node
net.neighbors <- get_neighbors_list(net)

# Number of model iterations to average the results (relates to ic_2_report)
n_iter <- 10


#### 3. Run the simulation, visualize results ####

ic_2_report <- IC_2_average(net.size, n_iter, net.neighbors)

# Prepare data to plot by averaging 
ic_2_plot <- ic_2_report %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.) %>%
  head(30)

# Plot average results:
data.frame(timestep=1:length(ic_2_plot), infected_percentage=ic_2_plot) %>%
  ggplot(., aes(timestep, infected_percentage)) + 
  geom_line() + geom_point(color="#56B4E9") + theme_minimal() +
  labs(
    title = "IC Algorithm (Case 2) (average of 100 iterations for all seeds)",
    x = "Day",
    y = "Fraction of the population who have seen the video"
  )


##### 3.1 Average for specific seed nodes ####

report_node_average <- pbapply::pblapply(
  1:net.size, 
  function(i) {
    replicate(100, IC_2(c(i), net.neighbors, net.size))
  }
)

spread_length <- colMeans(
  sapply(report_node_average, sapply, length)
)

res_df <- data.frame(
  node = 1:net.size,
  degree = degree(net),
  betwenness = betweenness(net),
  pshare = net.attributes$pYoutube,
  fully_spread = spread_length
) |>
  arrange(fully_spread)

View(rbind(head(res_df, 5),tail(res_df, 5)))


##### 3.2 Compare fastest and slowest seeds ####

set.seed(123)
IC_2_fastest <- replicate(
  100, IC_2(c(65), net.neighbors, net.size)
)

set.seed(123)
IC_2_slowest <- replicate(
  100, IC_2(c(78), net.neighbors, net.size)
)

ic_2_plot_fastest <- IC_2_fastest %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.)

ic_2_plot_slowest <- IC_2_slowest %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.)

list(fastest=ic_2_plot_fastest, slowest=ic_2_plot_slowest) %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  data.frame(
    timestep = 1:length(.[[1]])
  ) %>%
  pivot_longer(1:2, names_to = "speed", values_to = "contagion") %>%
  ggplot(aes(timestep, contagion, color=speed)) + geom_point() + 
  geom_line() + theme_minimal()
