# params_setting_v4.R

generate_JC_matrix <- function() {
  # Initialize matrix with 1/3 for all off-diagonal elements
  JC_matrix <- matrix(1/3, nrow = 5, ncol = 5)
  diag(JC_matrix) <- 0  # Diagonal elements are 0 for A, C, G, T
  JC_matrix[5, ] <- 0   # Transitions to 'N' are not allowed
  JC_matrix[, 5] <- 0   # Transitions from 'N' are not allowed
  JC_matrix[5, 5] <- 1  # 'N' transitions to itself with probability 1

  rownames(JC_matrix) <- colnames(JC_matrix) <- c("A", "C", "G", "T", "N")

  return(JC_matrix)
}

#' Generate Kimura 2-Parameter (K2P) Substitution Matrix
#'
#' @description
#' Generate a Kimura 2-Parameter (K2P) Transition Matrix
#'
#' This function generates a 5x5 transition probability matrix based on the Kimura 2-Parameter (K2P) model.
#' The matrix includes transition rates for nucleotides (A, C, G, T) and an additional state (N),
#' where transitions and transversions are governed by the parameters `alpha` and `beta`.
#'
#' @param alpha Numeric value representing the transition rate (rate of A<->G and C<->T mutations)
#' @param beta Numeric value representing the transversion rate (rate of all other mutations between A, C, G, T)
#'
#' @return A 5x5 matrix with row and column names (A, C, G, T, N) containing substitution
#' probabilities. The matrix entries represent the probability of transitioning from the row
#' nucleotide to the column nucleotide.
#'
#' @details
#' The K2P model assumes:
#' \itemize{
#'   \item Different rates for transitions (alpha) and transversions (beta)
#'   \item Equal base frequencies
#'   \item The 'N' state can only transition to itself
#'   \item No transitions are allowed between any nucleotide and 'N'
#' }
#' Matrix entries are normalized by (2*beta + alpha) to ensure proper probability distribution.
#'
#' @examples
#' # Generate K2P matrix with transition rate 0.3 and transversion rate 0.1
#' K2P <- generate_K2P_matrix(alpha = 0.3, beta = 0.1)
#'
#' @export
generate_K2P_matrix <- function(alpha, beta) {
  # Initialize a 5x5 matrix with specified transition and transversion rates
  K2P_matrix <- matrix(beta/(2*beta + alpha), nrow = 5, ncol = 5)

  # Set transition rates for A, C, G, T
  K2P_matrix[1, 3] <- alpha/(2*beta + alpha)  # A->G
  K2P_matrix[3, 1] <- alpha/(2*beta + alpha)  # G->A
  K2P_matrix[2, 4] <- alpha/(2*beta + alpha)  # C->T
  K2P_matrix[4, 2] <- alpha/(2*beta + alpha)  # T->C

  # Set diagonal to 0 for A, C, G, T
  diag(K2P_matrix) <- 0

  # Add the N state transitions
  K2P_matrix[5, ] <- c(0, 0, 0, 0, 1)  # N can only transition to N
  K2P_matrix[, 5] <- c(0, 0, 0, 0, 1)  # No other state can transition to N

  rownames(K2P_matrix) <- colnames(K2P_matrix) <- c("A", "C", "G", "T", "N")

  return(K2P_matrix)
}

