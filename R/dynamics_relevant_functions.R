# dynamics_relevant_functions_v5.R

sample_cells <- function(cell_info, num_samples) {
  # Filter cells that are not dead
  live_cells <- cell_info[is.na(cell_info$death_time), ]

  # If there are fewer live cells than the number of samples requested,
  # return a warning and sample all live cells
  if (nrow(live_cells) < num_samples) {
    warning("Number of samples requested is greater than the number of live cells. Sampling all live cells.")
    num_samples <- nrow(live_cells)
  }

  # Sample cells
  sampled_cells <- live_cells[sample(nrow(live_cells), num_samples), ]

  # Return cell_info of sampled cell
  return(sampled_cells)
}


sample_action <- function(total_birth_rate, total_death_rate, total_transition_rate) {
  # Combine rates into a vector
  rates <- c(total_birth_rate, total_death_rate, total_transition_rate)

  # Define action names
  actions <- c("birth", "death", "transition")

  # Sample an action based on its relative rate
  action <- sample(actions, size = 1, prob = rates)

  return(action)
}


pick_action_cell <- function(population, action_rate) {
  # Calculate the weighted probabilities
  weighted_rates <- population * action_rate

  # Check if all weighted rates are zero
  if (sum(weighted_rates) == 0) {
    stop("No available cells to pick from. All rates are zero.")
  }

  # Choose a clone proportional to its action rate
  cell <- sample(names(population), size = 1, prob = weighted_rates)

  return(cell)
}



pick_transition <- function(population, edge_transition_rates){
  # Choose a transition proportional to it's transition rate
  transition_parents <- edge_transition_rates$parent

  transition_idx <- sample(1:nrow(edge_transition_rates), size = 1, prob = population[transition_parents] * edge_transition_rates$rate)

  transition_clone <- edge_transition_rates[transition_idx, 'child']
  # update transtion matrix
  edge_transition_rates[transition_idx, 'rate'] = 0

  return(list(transition_clone = transition_clone, edge_transition_rates = edge_transition_rates))

}


plot_population_history_archived <- function(simulation_result, clone_colors) {
  require(ggplot2)
  require(reshape2)

  # Transform the data to long format
  long_data <- reshape2::melt(t(simulation_result$population_history))
  colnames(long_data) <- c("time", "clone", "population")

  # Create the plot
  p <- ggplot(long_data, aes(x = time, y = population, group  = clone)) +
    geom_line(aes(color = clone)) +
    labs(title = "Population History", x = "Time", y = "Population") +
    theme_bw() +
    scale_color_manual(values = clone_colors, name = "Clone") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title

  return(p)
}

plot_population_history <- function(simulation_result, clone_colors, title = "Population History") {
  require(ggplot2)
  require(reshape2)

  # Transform the data to long format
  long_data <- reshape2::melt(t(simulation_result$population_history))
  colnames(long_data) <- c("time", "clone", "population")

  p <- ggplot(long_data, aes(x = time, y = population, group = clone)) +
    geom_line(aes(color = clone)) +
    labs(title = title, x = "Time", y = "Population") +  # Use the title parameter here
    theme_bw() +
    scale_color_manual(values = clone_colors, name = "Clone") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),   # Centered, larger, and bold title
      axis.title.x = element_text(size = 14, face = "bold"),               # Bold x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),               # Bold y-axis title
      axis.text.x = element_text(size = 12),                               # Larger x-axis tick labels
      axis.text.y = element_text(size = 12),                               # Larger y-axis tick labels
      legend.title = element_text(size = 14, face = "bold"),               # Bold legend title
      legend.text = element_text(size = 12)                                # Larger legend text
    )


  return(p)
}

simulate_sc_dynamics <- function(initial_population, max_steps, intrinsic_birth_rates, intrinsic_death_rates, edge_transition_rates, clone_capacity, chr_lengths, mutation_rate) {
  time <- 0
  steps <- 0

  time_history <- time
  population <- initial_population
  population_history <- data.frame(t1 = population)

  # Initialize cell info data frame
  cell_info <- data.frame(clone = rep(names(initial_population), initial_population),
                          parent = NA,
                          birth_time = rep(0, sum(initial_population)),
                          death_time = NA,
                          cell_index = 1:sum(initial_population),
                          clone_index = as.numeric(unlist(lapply(initial_population, function(x) if (x > 0) seq_len(x) else integer(0)))),
                          stringsAsFactors = FALSE)


  mutation_info <- data.frame(clone = character(),
                              cell_index = integer(),
                              haplotype = character(),
                              chrom = character(),
                              pos = integer(),
                              time = numeric(),
                              stringsAsFactors = FALSE)

  while (steps < max_steps) {
    birth_rates <- pmax(intrinsic_birth_rates*(1 - sum(population)/clone_capacity), 0)  # Ensure rates are non-negative
    death_rates <- pmax(intrinsic_death_rates*population, 0)
    total_birth_rate <- sum(birth_rates)
    total_death_rate <- sum(death_rates)
    total_transition_rate <- sum(population[edge_transition_rates$parent]*edge_transition_rates$rate)
    total_rate <- total_birth_rate + total_death_rate + total_transition_rate

    # Check if total_rate is zero (i.e., all rates are zero)
    if (total_rate <= 0) {
      warning("Simulation stopped early due to zero total rate. Check your parameters.")
      break  # Exit the loop if no events can occur
    }

    # Step 3: Gillespie step
    time_step <- rexp(1, rate = total_rate)
    time <- time + time_step
    steps <- steps + 1

    action <- sample_action(total_birth_rate, total_death_rate, total_transition_rate)

    # Step 4: Update population
    if (action == 'birth') {
      birth_clone <- pick_action_cell(population, birth_rates)
      population[birth_clone] <- population[birth_clone] + 1
      # one cell to divide
      action_cell_idx <- select_cell_index(cell_info, birth_clone)
      # as cell divides, it will die
      cell_info$death_time[cell_info$cell_index == action_cell_idx] <- time

      # Record birth in cell info
      cell_index <- c(nrow(cell_info) + 1, nrow(cell_info) + 2)
      clone_index <- c(sum(cell_info$clone == birth_clone) + 1, sum(cell_info$clone == birth_clone) + 2)
      new_info <- data.frame(clone = birth_clone,
                             parent = action_cell_idx,
                             birth_time = time,
                             death_time = NA,
                             cell_index = cell_index,
                             clone_index = clone_index)

      cell_info <- rbind(cell_info, new_info)

      if (runif(1) < mutation_rate) {
        # Select a chromosome and position for the mutation
        clone_chr_lengths <- chr_lengths[[birth_clone]]
        mut_hap <- sample(names(clone_chr_lengths), 1)
        mut_chr <- sample(names(clone_chr_lengths[[mut_hap]]), 1, prob = clone_chr_lengths[[mut_hap]])
        mut_pos <- sample(1:clone_chr_lengths[[mut_hap]][mut_chr], 1)


        # Choose one of the daughter cells to inherit the mutation
        daughter_cell_idx <- sample(cell_index, 1)

        mutation_info <- rbind(mutation_info, data.frame(clone = birth_clone,
                                                         cell_index = daughter_cell_idx,
                                                         haplotype = mut_hap,
                                                         chrom = mut_chr,
                                                         pos = mut_pos,
                                                         time = time))

      }

    } else if (action == 'death') {
      dead_clone <- pick_action_cell(population, death_rates)
      action_cell_idx <- select_cell_index(cell_info, dead_clone)
      population[dead_clone] <- population[dead_clone] - 1
      # Record death in cell info
      cell_info$death_time[cell_info$cell_index == action_cell_idx] <- time

    } else if (action == 'transition') {
      transition_update <- pick_transition(population, edge_transition_rates)
      edge_transition_rates <- transition_update$edge_transition_rates # transition only happens once, the rate will turn to 0
      transition_clone <- transition_update$transition_clone

      population[transition_clone] <- population[transition_clone] + 1

      parent_clone = edge_transition_rates[edge_transition_rates$child == transition_clone,]$parent

      # one cell to transformation
      action_cell_idx <- select_cell_index(cell_info, parent_clone)

      # as cell transforms, it will die
      cell_info$death_time[cell_info$cell_index == action_cell_idx] <- time

      # Record transition in cell info
      cell_index <- c(nrow(cell_info) + 1, nrow(cell_info) + 2)
      clone_index <- c(sum(cell_info$clone == parent_clone) + 1, 1)

      new_info <- data.frame(clone = c(parent_clone, transition_clone),
                             parent = action_cell_idx,
                             birth_time = time,
                             death_time = NA,
                             cell_index = cell_index,
                             clone_index = clone_index)
      cell_info <- rbind(cell_info, new_info)


    }
    population_history <- cbind(population_history, population)
    time_history <- c(time_history, time)
  }
  colnames(population_history) <- names(time_history) <- paste0("t", 1:ncol(population_history))

  return(list(population_history = population_history, time_history = time_history, cell_info = cell_info, mutation_info = mutation_info))
}



select_cell_index <- function(cell_info, clone_type) {

  # Filter cell_info based on the clone type
  cell_indices <- cell_info[is.na(cell_info$death_time) & cell_info$clone == clone_type,]$cell_index

  # Select the cell index
  if (length(cell_indices) == 1) {
    return(cell_indices)
  } else {
    return(sample(cell_indices, 1))
  }
}
