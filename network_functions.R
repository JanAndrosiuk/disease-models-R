#### Functions required to run spread models ####

#' Creates a vector of neighbors for each node.
#' @param network network object of "igraph" class
#' @return list -> names as nodes id, values as vectors of neighbor ids.

get_neighbors_list <- function(network) {
  
  node_neighbors <- as.matrix(igraph::as_adjacency_matrix(network, type='both'))
  
  # Find connections (with reverse direction duplicates).
  node_neighbors <- which(node_neighbors>0, arr.ind=TRUE)
  
  # split second columns based on values of the first
  node_neighbors <- split(node_neighbors[,1], node_neighbors[,2])
  
  return(node_neighbors)
}


#' Performs non-greedy difference of two lists (sets).
#' @param set (set1, set2); vectors
#' @return vector
#' @example [1, 2, 2, 2, 3, 3], [2, 3] -> [1, 2, 2, 3]

setdiff_nongreedy <- function(set1, set2) {
  if (length(which(is.na(match(set1, set2)) == F)) != 0) {
    return(set1[-which(is.na(match(set1, set2)) == F)])
  } else {
    return(set1)
  }
}


#' Samples the elements of a vector based on a given probability.
#' @param vect target vector
#' @param proba uniform probability of each node getting sampled
#' @return sampled vector

sample_vect <- function(vect, proba) {
  return(
    vect[which(
      runif(length(vect), 0, 1) <= proba
    )]
  )
}

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
  
  # Return cumulated spread log.
  # return(cumsum(spread.log) / net.size)
  
  # Return last element of the cumulated sum. 
  return((cumsum(spread.log)/net.size) %>% .[length(.)])
}


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
  # return(cumsum(spread.log) / net.size)
  return(length(spread.log))
}
