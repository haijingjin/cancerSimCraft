% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamics_relevant_functions.R
\name{simulate_sc_dynamics}
\alias{simulate_sc_dynamics}
\title{Simulate Single-Cell Dynamics with Mutations}
\usage{
simulate_sc_dynamics(
  initial_population,
  max_steps,
  intrinsic_birth_rates,
  intrinsic_death_rates,
  edge_transition_rates,
  clone_capacity,
  chr_lengths,
  mutation_rate
)
}
\arguments{
\item{initial_population}{Named numeric vector specifying the initial count for each clone type}

\item{max_steps}{Maximum number of simulation steps to perform}

\item{intrinsic_birth_rates}{Named numeric vector of birth rates for each clone type}

\item{intrinsic_death_rates}{Named numeric vector of death rates for each clone type}

\item{edge_transition_rates}{Data frame with columns: parent, child, rate - specifying transitions between clone types}

\item{clone_capacity}{Numeric value representing the carrying capacity of the system of each clone}

\item{chr_lengths}{Nested list structure defining chromosome lengths for each clone and haplotype}

\item{mutation_rate}{Probability of mutation occurring during cell division}
}
\value{
A list containing:
\itemize{
\item population_history: Data frame tracking population counts over time
\item time_history: Vector of time points corresponding to population_history
\item cell_info: Data frame with detailed lineage information for each cell
\item mutation_info: Data frame recording mutation events that occurred during simulation
}
}
\description{
Simulates cell population dynamics using a Gillespie algorithm, tracking births, deaths,
transitions between cell types, and mutations on chromosomes. This function models cell
growth with logistic constraints and generates detailed lineage information.
}
\details{
The simulation implements a Gillespie algorithm with three possible events:
\enumerate{
\item Birth: A cell divides into two daughter cells with possible mutations
\item Death: A cell dies and is removed from the population
\item Transition: A cell changes from one clone type to another
}

Birth rates are modulated by logistic growth constraints based on total population
and clone capacity. The simulation tracks detailed lineage information including parent-child
relationships, birth and death times for each cell, and mutation events.
}
\note{
The simulation stops when either max_steps is reached or the total event rate becomes zero.
}
\seealso{
\code{\link{sample_action}}, \code{\link{pick_action_cell}}, \code{\link{select_cell_index}}, \code{\link{pick_transition}}
}
