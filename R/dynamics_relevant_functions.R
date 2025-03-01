#' Sample Living Cells from Simulation Results
#'
#' @description
#' Randomly samples a specified number of living cells from the cell lineage information
#' produced by a dynamics simulation.
#'
#' @param cell_info Data frame containing cell lineage information with columns:
#'        clone, parent, birth_time, death_time, cell_index, clone_index
#' @param num_samples Integer specifying the number of cells to sample
#'
#' @return A data frame containing the cell_info rows for the randomly sampled cells
#'
#' @details
#' This function first filters the cell_info data frame to identify living cells
#' (those with NA in death_time), then randomly samples the specified number of cells
#' from this subset.
#'
#' If the number of requested samples exceeds the number of living cells, the function
#' issues a warning and returns all living cells instead.
#'
#' @examples
#' # Example usage:
#' cell_info <- data.frame(
#'   clone = c("Clone1", "Clone1", "Clone2"),
#'   death_time = c(NA, 5, NA),
#'   cell_index = c(1, 2, 3)
#' )
#' num_samples <- 2
#'
#' sampled_cells <- sample_cells(cell_info, num_samples)
#' print(sampled_cells)  # Sampled live cells
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}
#'
#' @export
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

#' Sample a Gillespie Event Type
#'
#' @description
#' Randomly selects an event type (birth, death, or transition) based on the relative rates
#' of each event in a Gillespie simulation algorithm.
#'
#' @param total_birth_rate Numeric value representing the sum of all birth rates in the system
#' @param total_death_rate Numeric value representing the sum of all death rates in the system
#' @param total_transition_rate Numeric value representing the sum of all transition rates in the system
#'
#' @return A character string: either "birth", "death", or "transition" indicating the selected event
#'
#' @details
#' This function implements the event selection step of the Gillespie algorithm.
#' The probability of selecting each event type is proportional to its rate
#' relative to the total rate of all events combined.
#'
#' @examples
#' # Example usage:
#' total_birth_rate <- 0.5
#' total_death_rate <- 0.3
#' total_transition_rate <- 0.2
#'
#' action <- sample_action(total_birth_rate, total_death_rate, total_transition_rate)
#' print(action)  # Possible output: "birth"
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}
#'
#'
#' @export
sample_action <- function(total_birth_rate, total_death_rate, total_transition_rate) {
  # Combine rates into a vector
  rates <- c(total_birth_rate, total_death_rate, total_transition_rate)
  # Define action names
  actions <- c("birth", "death", "transition")
  # Sample an action based on its relative rate
  action <- sample(actions, size = 1, prob = rates)
  return(action)
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

#' Select a Clone Type for an Event
#'
#' @description
#' Selects a clone type on which an event (birth, death, etc.) will occur based
#' on the population size and action rate of each clone type.
#'
#' @param population Named numeric vector containing the current population count for each clone type
#' @param action_rate Named numeric vector containing the per-cell rate for the action for each clone type
#'
#' @return A character string representing the name of the selected clone type
#'
#' @details
#' This function implements the clone type selection step in the Gillespie algorithm.
#' The probability of selecting a particular clone type is proportional to its
#' population multiplied by its action rate. This ensures that more abundant clones
#' and clones with higher rates are more likely to be selected.
#'
#' @examples
#' # Example usage:
#' population <- c(Clone1 = 10, Clone2 = 5)
#' action_rate <- c(Clone1 = 0.1, Clone2 = 0.05)
#'
#' selected_clone <- pick_action_cell(population, action_rate)
#' print(selected_clone)  # Possible output: "Clone1"
#'
#' @seealso
#' \code{\link{sample_action}}, \code{\link{simulate_sc_dynamics}}
#'
#' @export
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


#' Select a Transition Event Between Clone Types
#'
#' @description
#' Selects a specific transition event between clone types based on the current
#' population and edge transition rates, then updates the transition rate matrix
#' to prevent the same transition from occurring again.
#'
#' @param population Named numeric vector containing the current population count for each clone type
#' @param edge_transition_rates Data frame with columns: parent, child, rate - specifying transition probabilities
#'        between clone types
#'
#' @return A list containing:
#' \itemize{
#'   \item transition_clone: Character string representing the name of the target clone type (the clone to which the cell transitions)
#'   \item edge_transition_rates: Updated data frame of transition rates with the selected
#'         transition's rate set to zero
#' }
#'
#' @details
#' This function implements the transition selection step in the Gillespie algorithm for
#' cell clone type transitions. The probability of selecting a particular transition is
#' proportional to the population of the parent clone multiplied by the transition rate.
#' After selection, the transition rate for the selected transition is set to zero,
#' indicating it can only occur once in the simulation.
#'
#' @examples
#' # Population counts
#' pop <- c(A = 100, B = 50, C = 30)
#'
#' # Define transition matrix
#' transitions <- data.frame(
#'   parent = c("A", "B", "A"),
#'   child = c("B", "C", "C"),
#'   rate = c(0.01, 0.02, 0.005)
#' )
#'
#' # Select and process a transition event
#' result <- pick_transition(pop, transitions)
#' result$transition_clone  # The selected target clone
#' result$edge_transition_rates  # Updated transition matrix
#'
#' @seealso
#' \code{\link{sample_action}}, \code{\link{pick_action_cell}}, \code{\link{simulate_sc_dynamics}}
#'
#'
#' @export
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

#' Plot Population Dynamics from Simulation Results
#'
#' @description
#' Creates a ggplot visualization of clone population dynamics over time from
#' simulation results produced by the simulate_sc_dynamics function.
#'
#' @param simulation_result List object returned by simulate_sc_dynamics containing population_history
#' @param clone_colors Named vector of colors where names match clone types in the simulation
#' @param title Character string for the plot title, defaults to "Population History"
#'
#' @return A ggplot object displaying population dynamics over time with customized styling
#'
#' @details
#' This function transforms the population history data from wide to long format and
#' creates a line plot showing how the population of each clone type changes over time.
#' The plot includes custom styling for readability and aesthetic presentation.
#'
#' The time points displayed on the x-axis correspond to the simulation steps without
#' displaying the raw time values (using discrete scale with no labels).
#'
#' @examples
#'
#' # Define colors for each clone type
#' clone_colors <- c(A = "blue", B = "red")
#'
#' # Create and display the plot
#' p <- plot_population_history(sim_result, clone_colors, "Clone Evolution")
#' print(p)
#'
#' # Save the plot
#' ggsave("population_dynamics.png", p, width = 8, height = 6)
#'
#' @seealso
#' \code{\link{simulate_sc_dynamics}}
#'
#' @import ggplot2
#' @import reshape2
#'
#' @export
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

#' Simulate Single-Cell Dynamics with Mutations
#'
#' @description
#' Simulates cell population dynamics using a Gillespie algorithm, tracking births, deaths,
#' transitions between cell types, and mutations on chromosomes. This function models cell
#' growth with logistic constraints and generates detailed lineage information.
#'
#' @param initial_population Named numeric vector specifying the initial count for each clone type
#' @param max_steps Maximum number of simulation steps to perform
#' @param intrinsic_birth_rates Named numeric vector of birth rates for each clone type
#' @param intrinsic_death_rates Named numeric vector of death rates for each clone type
#' @param edge_transition_rates Data frame with columns: parent, child, rate - specifying transitions between clone types
#' @param clone_capacity Numeric value representing the carrying capacity of the system of each clone
#' @param chr_lengths Nested list structure defining chromosome lengths for each clone and haplotype
#' @param mutation_rate Probability of mutation occurring during cell division
#'
#' @return A list containing:
#' \itemize{
#'   \item population_history: Data frame tracking population counts over time
#'   \item time_history: Vector of time points corresponding to population_history
#'   \item cell_info: Data frame with detailed lineage information for each cell
#'   \item mutation_info: Data frame recording mutation events that occurred during simulation
#' }
#'
#' @details
#' The simulation implements a Gillespie algorithm with three possible events:
#' 1. Birth: A cell divides into two daughter cells with possible mutations
#' 2. Death: A cell dies and is removed from the population
#' 3. Transition: A cell changes from one clone type to another
#'
#' Birth rates are modulated by logistic growth constraints based on total population
#' and clone capacity. The simulation tracks detailed lineage information including parent-child
#' relationships, birth and death times for each cell, and mutation events.
#'
#' @note
#' The simulation stops when either max_steps is reached or the total event rate becomes zero.
#'
#'
#' @seealso
#' \code{\link{sample_action}}, \code{\link{pick_action_cell}}, \code{\link{select_cell_index}}, \code{\link{pick_transition}}
#'
#'
#' @export
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


#' Select a Specific Cell from a Clone Population
#'
#' @description
#' Randomly selects a specific cell from the living cells of a given clone type
#' based on the cell lineage information.
#'
#' @param cell_info Data frame containing cell lineage information with columns:
#'        clone, parent, birth_time, death_time, cell_index, clone_index
#' @param clone_type Character string specifying the clone type from which to select a cell
#'
#' @return An integer representing the selected cell's unique index
#'
#' @details
#' This function filters the cell_info data frame to find all living cells (those without a death_time)
#' of the specified clone type. If there is only one such cell, it is returned directly.
#' If there are multiple cells, one is chosen randomly with equal probability.
#'
#' Living cells are identified as those where death_time is NA in the cell_info data frame.
#'
#' # Example usage:
#' cell_info <- data.frame(
#'   clone = c("Clone1", "Clone1", "Clone2"),
#'   death_time = c(NA, 5, NA),
#'   cell_index = c(1, 2, 3)
#' )
#' clone_type <- "Clone1"
#'
#' selected_cell_index <- select_cell_index(cell_info, clone_type)
#' print(selected_cell_index)  # Possible output: 1
#'
#' @seealso
#' \code{\link{pick_action_cell}}, \code{\link{simulate_sc_dynamics}}
#'
#'
#' @export
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
