% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamics_relevant_functions.R
\name{select_cell_index}
\alias{select_cell_index}
\title{Select a Specific Cell from a Clone Population}
\usage{
select_cell_index(cell_info, clone_type)
}
\arguments{
\item{cell_info}{Data frame containing cell lineage information with columns:
clone, parent, birth_time, death_time, cell_index, clone_index}

\item{clone_type}{Character string specifying the clone type from which to select a cell}
}
\value{
An integer representing the selected cell's unique index
}
\description{
Randomly selects a specific cell from the living cells of a given clone type
based on the cell lineage information.
}
\details{
This function filters the cell_info data frame to find all living cells (those without a death_time)
of the specified clone type. If there is only one such cell, it is returned directly.
If there are multiple cells, one is chosen randomly with equal probability.

Living cells are identified as those where death_time is NA in the cell_info data frame.
}
\section{Example usage:}{
cell_info <- data.frame(
clone = c("Clone1", "Clone1", "Clone2"),
death_time = c(NA, 5, NA),
cell_index = c(1, 2, 3)
)
clone_type <- "Clone1"

selected_cell_index <- select_cell_index(cell_info, clone_type)
print(selected_cell_index)  # Possible output: 1
}

\seealso{
\code{\link{pick_action_cell}}, \code{\link{simulate_sc_dynamics}}
}
