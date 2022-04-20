#### Threshold Model Algorithm ####

# ASSUMPTION 1:
  # This model is deterministic, not stochastic

# ASSUMPTION 2:
  # Nodes become active when the following condition becomes true:
  # number_of_neighbors * node_i_fractional_thres < number_of_active_neighbors

# ASSUMPTION 3:
  # Based on #A2 -> this model won't necessarily reach 100% contagion for 
  # each set of seed nodes.

#### 0. Libraries, functions ####

rm(list=ls())

library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path)) 

source("network_functions.R")

# library(visNetwork)
# library(htmlwidgets)
library(igraph)
library(dplyr)
library(zoo) # locf function
library(tidyverse)
library(pbapply)


#### 1. Simulation  function ####

threshold_model <- function(spread.seeds, net.neighbors, net.size){
  # Save number of newly activated nodes for each day in the log.
  spread.log <- c()
  
  # Column of the network attributes containing fractional threshold.
  thres_attr <- "tRecipe"
  
  # Day 0
  
  # spread.status informs about nodes that were activated during the simulation.
  # Initially every nodes is inactive (State=0).
  spread.status = rep(0, net.size)
  
  # Day 1
  
  spread.status[spread.seeds] <- 1
  spread.log[[1]] <- length(spread.seeds)
  
  # Day 2+
  stop_criterion <- T
  while (stop_criterion) {
    
    # Get list of active nodes which will have a chance to spread phenomena
    # across their neighbors.
    spread.active = which(spread.status == 1)
    
    # 1) Find nodes nodes that are unactive, but have at least one active neighbor.
    # 2) Save all neighbors of those nodes. 
    # 3) Save active neighbors of those nodes.
    # 4) Get thresholds of target nodes
    spread.targets <- setdiff(
      unlist(net.neighbors[spread.active], use.names=F),
      which(spread.status == 1)
    ) 
    
    spread.all_neighbors <- net.neighbors[spread.targets]
    
    spread.active_neighbors <- unlist(lapply(
      lapply(
        spread.all_neighbors, setdiff, which(spread.status == 0)
      ), length
    ), use.names = F)
    
    spread.all_neighbors <- unlist(lapply(
      spread.all_neighbors, length
    ), use.names = F)
    
    spread.thresholds <- net.attributes[,thres_attr][spread.targets]
    
    # Nodes are infected given the condition: N_i * threshold_i < N_active_i
    # where N_i and N_active_i stand for all neighbors and active neighbors
    # of target nodes.
    spread.new <- spread.targets[
      spread.all_neighbors * spread.thresholds < spread.active_neighbors
    ]
    
    # If there are any new infected cases, change their status, and save their
    # occurrence to the log.
    if (length(spread.new) > 0) {
      spread.status[spread.new] <- 1
      spread.log[[length(spread.log)+1]] <- length(spread.new)  
    } else {
      spread.log[[length(spread.log)+1]] <- 0
    }
    # print(sum(spread.status)/net.size)
    stop_criterion <- sum(unlist(tail(spread.log, 5))) != 0
  }
  # print(spread.log)
  # print(unlist(spread.log, use.names = F))
  return(cumsum(spread.log) / net.size)
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

# Define fractional threshold column name
thres_attr = "tRecipe"

#### 3. Run the simulation, visualize results ####

threshold_model(4, net.neighbors, net.size)

thres_report <- pbapply::pblapply(
  1:net.size, threshold_model, net.neighbors, net.size
)

##### 3.1 Plot average path for different seeds ####

# Prepare data to plot by averaging 
thres_plot <- thres_report %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.) %>%
  head(30)

# Plot average results:
data.frame(timestep=1:length(thres_plot), infected_percentage=thres_plot) %>%
  ggplot(., aes(timestep, infected_percentage)) + 
  geom_line() + geom_point(color="#56B4E9") + theme_minimal() +
  labs(
    title = "Threshold model (average of all seeds)",
    x = "Day",
    y = "Fraction of the population who have seen the video"
  )


##### 3.2 Check best seeds ####
report_node_average <- pbapply::pblapply(
  1:net.size, threshold_model, net.neighbors, net.size
) 

spread_max <- sapply(report_node_average, max)

res_df <- data.frame(
  node = 1:net.size,
  degree = degree(net),
  betwenness = betweenness(net),
  pshare = net.attributes[,thres_attr],
  max_spread = spread_max
) |>
  arrange(desc(max_spread))

View(rbind(head(res_df, 5),tail(res_df, 5)))
