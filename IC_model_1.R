#### Independent Cascade Model with internally equal probabilities ####

# ASSUMPTION 1:
  # It doesn't matter if we analyze infectious nodes one by one
  # as the probability of neighbors getting infected is constant.

# ASSUMPTION 2:
  # Infectious are only infectious only for specific range of days.
  # Refers to the infectious period.

# ASSUMPTION 3:
  # Algorithm will run until the total number of infected (sum(inf.status))
  # won't change for number of time steps equal to the infectious period.


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


#### 1. Simulation function ####
IC <- function(inf.seed, inf.proba, inf.inf_period, net.neighbors) {
  
  # Save infection log for each day
  inf.log <- c()
  
  # Day 0
  inf.status = rep(0, net.size)
  
  # Day 1
  inf.status[inf.seed] <- 1
  inf.log[[length(inf.log)+1]] <- inf.seed
  
  # Day 2 till the end
  
  # stop_criterion refers to the initial value of the argument which is used to 
  # perform the while loop. Check description at the end of the loop.
  stop_criterion = T
  
  while (stop_criterion) {
    
    # get infected nodes which are still infectious
    # based on the infectious period
    inf.current_infectious <- unlist(tail(inf.log, inf.inf_period))
    
    # get neighbors of currently infected
    inf.current_neighbors <- setdiff(
      unlist(net.neighbors[inf.current_infectious], use.names=F),
      which(inf.status == 1)
    )
    
    # from the neighbors above, get those who got infected
    inf.infected_neighbors <- sample_vect(
      inf.current_neighbors, inf.proba
    )
    
    if (length(inf.infected_neighbors) > 0) {
      
      # Updating infection status of population
      inf.status[inf.infected_neighbors] <- 1
      
      # save currently infected neighbors to infection log
      # they will be infectious from in the next 3 days
      inf.log[[length(inf.log)+1]] <- inf.infected_neighbors  
    }
    else {
      inf.log[[length(inf.log)+1]] <- integer(0)
    }
    
    # check if number of infected nodes (since beginning of simulation)
    # has changed over last n iterations (with n set to inf_period).
    # If it changed, continue simulation, if hasn't - stop.
    # TRUE indicates that it had changed, FALSE that didn't
    stop_criterion <- sum(sapply(tail(inf.log, inf.inf_period), length)) != 0
  }
  
  return(cumsum(sapply(inf.log, length)/140))
}


IC_average <- function(net.size, no_reps, inection_proba, infection_period,
                       network_neighbors) {
  return(
    pbapply::pbreplicate(
      no_reps,
      lapply(
        1:net.size, IC, inection_proba, infection_period, network_neighbors
      )
    )
  )
}


#### 2. Simulation parameters ####

# Load network object
net <- readRDS("data/net1.rds")
net.size = igraph::gorder(net)

# Define infection seed statically
inf.seed = 140

# Infection probability
inf.proba = 0.1

# Set the period of being infectious after being infected
inf.inf_period = 3

# Get neighbor vector for each node
net.neighbors <- get_neighbors_list(net)

# Number of model iterations to average the results (relates to ic_report)
n_iter <- 10


#### 3. Run the simulation, visualize results ####

ic_report <- IC_average(
  net.size, n_iter, inf.proba, inf.inf_period, net.neighbors
)

ic_plot <- ic_report %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.)
  
data.frame(
  timestep=1:length(ic_plot), 
  infected_percentage=ic_plot
) %>%
  ggplot(., aes(timestep, infected_percentage)) + 
    geom_line() + geom_point(color="#56B4E9") + theme_minimal() +
    labs(
      title = "IC Algorithm (average of 100 iterations for all seeds)",
      x = "Day",
      y = "Infected fraction of the population"
  )


#### 4. Delete nodes with highest and lowest betweenness metric ####

net_del_high <- net %>% delete_vertices(., tail(order(betweenness(.)), 3))
net_del_low <- net %>% delete_vertices(., head(order(betweenness(.)), 3))

net_del_high.neighbors <- get_neighbors_list(net_del_high)
net_del_low.neighbors <- get_neighbors_list(net_del_low)

report_del_high <- IC_average(
  gorder(net_del_high), 100, inf.proba, inf.inf_period, net_del_high.neighbors
)

report_del_low <- IC_average(
  gorder(net_del_low), 100, inf.proba, inf.inf_period, net_del_low.neighbors
)

ic_plot_high <- report_del_high %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.)

ic_plot_low <- report_del_low %>%
  lapply(., function(x) {
    length(x) <- max(lengths(.)); x 
  }) %>%
  lapply(., na.locf) %>%
  do.call(rbind, .) %>%
  colMeans(.)

data.frame(
  timestep=1:length(ic_plot_high), 
  infected_percentage=ic_plot_high
) %>%
  ggplot(., aes(
    timestep, infected_percentage
  )) + 
    geom_line() + geom_point(color="#56B4E9") + theme_minimal() +
    labs(
      title = "IC Algorithm without 3 nodes of the highest betweenness",
      x = "Day",
      y = "Infected fraction of the population"
    )

data.frame(
  timestep=1:length(ic_plot_low), 
  infected_percentage=ic_plot_low
) %>%
  ggplot(., aes(
    timestep, infected_percentage
  )) + 
    geom_line() + geom_point(color="#56B4E9") + theme_minimal() +
    labs(
      title = "IC Algorithm without 3 nodes of the lowest betweenness",
      x = "Day",
      y = "Infected fraction of the population"
    )

##### 4.1 Combine plots ####

lowest_df_size <- min(c(
  length(ic_plot), length(ic_plot_high), length(ic_plot_low)
))

res_df <- data.frame(
  day = 1:lowest_df_size,
  all_vertices = head(ic_plot, lowest_df_size),
  without_highest = head(ic_plot_high, lowest_df_size),
  without_lowest = head(ic_plot_low, lowest_df_size)
) |>
  head(40) |>
  pivot_longer(cols=2:4, names_to = "vertices", values_to = "values")


ggplot(res_df, aes(day, values, color=vertices)) + 
  geom_line() + geom_point() + theme_minimal() +
  labs(
    title = "IC Algorithm for 3 different sets of nodes",
    x = "Day",
    y = "Infected fraction of the population"
  )
